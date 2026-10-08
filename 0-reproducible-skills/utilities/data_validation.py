"""
Common data validation and manifest utilities for reproducible pipelines.

Provides:
  - SHA-256 verification for downloaded files
  - Data manifest generation and validation
  - Shape/encoding checks for CSVs and arrays
  - Logging of data parsing decisions (what was excluded, why)
"""

import hashlib
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Union
import csv


log = logging.getLogger(__name__)


class DataManifest:
    """Record and validate file checksums for reproducibility."""
    
    def __init__(self, manifest_path: Union[str, Path]):
        """
        Initialize manifest.
        
        Parameters
        ----------
        manifest_path : str or Path
            Path to manifest JSON file.
        """
        self.path = Path(manifest_path)
        self.data = {}
        
        if self.path.exists():
            with open(self.path) as f:
                self.data = json.load(f)
    
    def add(self, file_path: Union[str, Path], description: str = "", compute: bool = True):
        """
        Add or update a file in the manifest.
        
        Parameters
        ----------
        file_path : str or Path
            Path to the file.
        description : str
            Human-readable description of the file.
        compute : bool
            If True, compute SHA-256 of the file. If False, just record the path.
        
        Returns
        -------
        str
            The SHA-256 hash (or empty string if compute=False).
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            log.warning(f"File does not exist: {file_path}")
            return ""
        
        sha = ""
        if compute:
            sha = self._compute_sha256(file_path)
        
        self.data[str(file_path)] = {
            "description": description,
            "size_bytes": file_path.stat().st_size,
            "sha256": sha,
        }
        
        return sha
    
    def validate(self, file_path: Union[str, Path]) -> bool:
        """
        Check that a file matches its recorded checksum.
        
        Parameters
        ----------
        file_path : str or Path
            Path to the file.
        
        Returns
        -------
        bool
            True if checksum matches or file has no recorded checksum.
        
        Raises
        ------
        ValueError
            If file is not in manifest or checksum mismatches.
        """
        file_path = Path(file_path)
        key = str(file_path)
        
        if key not in self.data:
            raise ValueError(f"File not in manifest: {file_path}")
        
        recorded = self.data[key].get("sha256", "")
        if not recorded:
            log.warning(f"No checksum recorded for {file_path}; skipping validation")
            return True
        
        actual = self._compute_sha256(file_path)
        if actual != recorded:
            raise ValueError(
                f"Checksum mismatch for {file_path}\n"
                f"  Expected: {recorded}\n"
                f"  Got:      {actual}"
            )
        
        log.info(f"Checksum validated: {file_path}")
        return True
    
    def save(self):
        """Write manifest to disk."""
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.path, 'w') as f:
            json.dump(self.data, f, indent=2)
        log.info(f"Manifest saved: {self.path}")
    
    @staticmethod
    def _compute_sha256(file_path: Path) -> str:
        """Compute SHA-256 hash of a file."""
        sha = hashlib.sha256()
        with open(file_path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b''):
                sha.update(chunk)
        return sha.hexdigest()


class CSVValidator:
    """Validate CSV files during parsing."""
    
    @staticmethod
    def check_encoding(file_path: Union[str, Path], expected_encodings: List[str] = None):
        """
        Detect the encoding of a CSV file.
        
        Parameters
        ----------
        file_path : str or Path
            Path to the file.
        expected_encodings : list, optional
            List of expected encodings (e.g., ['utf-8', 'latin-1']).
        
        Returns
        -------
        str
            The detected encoding.
        
        Raises
        ------
        ValueError
            If detected encoding is not in expected_encodings (if provided).
        """
        if expected_encodings is None:
            expected_encodings = ['utf-8', 'utf-8-sig', 'latin-1', 'iso-8859-1']
        
        for encoding in expected_encodings:
            try:
                with open(file_path, encoding=encoding) as f:
                    f.read(1024)
                log.info(f"File encoding: {encoding}")
                return encoding
            except (UnicodeDecodeError, LookupError):
                continue
        
        raise ValueError(f"Could not detect encoding for {file_path}")
    
    @staticmethod
    def check_shape(
        file_path: Union[str, Path],
        expected_columns: Optional[List[str]] = None,
        expected_rows: Optional[int] = None,
        encoding: str = 'utf-8',
    ):
        """
        Check the shape and column names of a CSV file.
        
        Parameters
        ----------
        file_path : str or Path
            Path to the CSV file.
        expected_columns : list, optional
            Expected column names. If provided, will be checked.
        expected_rows : int, optional
            Expected number of rows (excluding header).
        encoding : str
            File encoding.
        
        Returns
        -------
        dict
            Shape information: {rows, columns, missing_columns}.
        
        Raises
        ------
        ValueError
            If shape doesn't match expectations.
        """
        with open(file_path, encoding=encoding) as f:
            reader = csv.DictReader(f)
            actual_cols = reader.fieldnames or []
            row_count = sum(1 for _ in reader)
        
        info = {
            'rows': row_count,
            'columns': len(actual_cols),
            'column_names': actual_cols,
            'missing_columns': [],
        }
        
        if expected_columns:
            missing = set(expected_columns) - set(actual_cols)
            if missing:
                info['missing_columns'] = sorted(missing)
                raise ValueError(
                    f"Missing columns in {file_path}: {missing}"
                )
        
        if expected_rows is not None and row_count != expected_rows:
            raise ValueError(
                f"Row count mismatch in {file_path}: "
                f"expected {expected_rows}, got {row_count}"
            )
        
        return info


class ParseLog:
    """Log data parsing decisions for transparency and reproducibility."""
    
    def __init__(self, log_path: Union[str, Path]):
        """
        Initialize parse log.
        
        Parameters
        ----------
        log_path : str or Path
            Where to save the log.
        """
        self.log_path = Path(log_path)
        self.entries = []
    
    def add(self, stage: str, key: str, value: str, detail: str = ""):
        """
        Log a parsing decision.
        
        Parameters
        ----------
        stage : str
            Analysis stage (e.g., 'ingest', 'preprocess').
        key : str
            What was decided (e.g., 'genes_dropped', 'cells_filtered').
        value : str
            Quantitative result (e.g., '342', '89%').
        detail : str
            Human-readable explanation.
        """
        entry = {
            'stage': stage,
            'key': key,
            'value': value,
            'detail': detail,
        }
        self.entries.append(entry)
        log.info(f"[{stage}] {key}: {value} {detail}")
    
    def save(self):
        """Write parse log to disk as JSON."""
        self.log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.log_path, 'w') as f:
            json.dump(self.entries, f, indent=2)
        log.info(f"Parse log saved: {self.log_path}")
