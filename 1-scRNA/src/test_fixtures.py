"""
Reusable test fixtures and validation functions for reproducible pipelines.

Provides:
  - Fixture generators for common data types (YAML configs, mock matrices, etc.)
  - Validation helpers (data shape/dtype checks, NaN/inf detection)
  - Comparison utilities (numerical tolerance, set equality)

These are domain-agnostic and can be imported into any project's test suite.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Union, Tuple


class MockData:
    """Factory for generating test data without external dependencies."""
    
    @staticmethod
    def sparse_matrix(
        nrows: int,
        ncols: int,
        density: float = 0.1,
        dtype: np.dtype = np.float32,
        seed: int = 42,
    ):
        """
        Generate a random sparse matrix (scipy.sparse.csr_matrix).
        
        Parameters
        ----------
        nrows, ncols : int
            Matrix dimensions.
        density : float
            Fraction of non-zero entries (0-1).
        dtype : dtype
            Data type.
        seed : int
            Random seed.
        
        Returns
        -------
        scipy.sparse.csr_matrix
        """
        from scipy.sparse import random as sparse_random
        rng = np.random.RandomState(seed)
        return sparse_random(nrows, ncols, density=density, dtype=dtype, random_state=rng)
    
    @staticmethod
    def dataframe(
        nrows: int,
        ncols: int = 5,
        column_types: dict = None,
        seed: int = 42,
    ) -> pd.DataFrame:
        """
        Generate a random DataFrame.
        
        Parameters
        ----------
        nrows : int
            Number of rows.
        ncols : int
            Number of columns (used only if column_types is None).
        column_types : dict, optional
            Maps column name to generator function. E.g.:
            {'gene': lambda n: [f'GENE_{i}' for i in range(n)],
             'count': lambda n: np.random.randint(0, 100, n)}
        seed : int
            Random seed.
        
        Returns
        -------
        pd.DataFrame
        """
        rng = np.random.RandomState(seed)
        
        if column_types is None:
            data = {f'col_{i}': rng.randn(nrows) for i in range(ncols)}
        else:
            data = {name: gen(nrows) for name, gen in column_types.items()}
        
        return pd.DataFrame(data)
    
    @staticmethod
    def config_yaml(path: Union[str, Path], **kwargs):
        """
        Write a minimal YAML config file.
        
        Parameters
        ----------
        path : str or Path
            Where to write the file.
        **kwargs
            Config keys and values. Defaults to a minimal valid config.
        
        Returns
        -------
        Path
            Path to the written file.
        """
        import yaml
        
        defaults = {
            'seed': 42,
            'paths': {
                'raw': 'data/raw',
                'processed': 'data/processed',
                'results': 'results',
            },
            'stages': {},
        }
        defaults.update(kwargs)
        
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(path, 'w') as f:
            yaml.dump(defaults, f)
        
        return path


class Validator:
    """Validation helpers for test assertions."""
    
    @staticmethod
    def assert_valid_matrix(X, min_rows: int = 0, max_rows: int = None,
                            min_cols: int = 0, max_cols: int = None,
                            dtype=None):
        """
        Check that a matrix has valid shape and dtype.
        
        Parameters
        ----------
        X : array-like or sparse
            Matrix to validate.
        min_rows, max_rows, min_cols, max_cols : int
            Expected shape bounds (None means unchecked).
        dtype : dtype, optional
            Expected dtype.
        
        Raises
        ------
        AssertionError
            If validation fails.
        """
        from scipy.sparse import spmatrix
        
        is_sparse = isinstance(X, spmatrix)
        nrows, ncols = X.shape
        
        if min_rows is not None:
            assert nrows >= min_rows, f"Expected >= {min_rows} rows, got {nrows}"
        if max_rows is not None:
            assert nrows <= max_rows, f"Expected <= {max_rows} rows, got {nrows}"
        if min_cols is not None:
            assert ncols >= min_cols, f"Expected >= {min_cols} cols, got {ncols}"
        if max_cols is not None:
            assert ncols <= max_cols, f"Expected <= {max_cols} cols, got {ncols}"
        
        if dtype is not None and not is_sparse:
            assert X.dtype == dtype, f"Expected dtype {dtype}, got {X.dtype}"
    
    @staticmethod
    def assert_no_missing(X: Union[np.ndarray, pd.DataFrame, pd.Series]):
        """
        Check for NaN, inf, or missing values.
        
        Parameters
        ----------
        X : array-like or DataFrame
            Data to check.
        
        Raises
        ------
        AssertionError
            If any missing values found.
        """
        if isinstance(X, (pd.DataFrame, pd.Series)):
            missing = X.isna().sum().sum() if isinstance(X, pd.DataFrame) else X.isna().sum()
            assert missing == 0, f"Found {missing} missing values"
        else:
            X = np.asarray(X)
            assert not np.isnan(X).any(), "Found NaN values"
            assert not np.isinf(X).any(), "Found inf values"
    
    @staticmethod
    def assert_close(a: np.ndarray, b: np.ndarray, rtol: float = 1e-5, atol: float = 1e-8):
        """
        Check that two arrays are numerically close.
        
        Parameters
        ----------
        a, b : array-like
            Arrays to compare.
        rtol, atol : float
            Relative and absolute tolerance.
        
        Raises
        ------
        AssertionError
            If arrays are not close.
        """
        a, b = np.asarray(a), np.asarray(b)
        assert a.shape == b.shape, f"Shape mismatch: {a.shape} vs {b.shape}"
        
        close = np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True)
        if not close:
            rel_err = np.abs((a - b) / (np.abs(b) + atol))
            max_err = np.nanmax(rel_err)
            raise AssertionError(f"Arrays not close; max relative error: {max_err}")


class Compare:
    """Comparison utilities for test assertions."""
    
    @staticmethod
    def sets_equal(a: set, b: set, name: str = "set"):
        """
        Assert two sets are equal, with detailed diff on failure.
        
        Parameters
        ----------
        a, b : set
            Sets to compare.
        name : str
            Name for error message.
        
        Raises
        ------
        AssertionError
            If sets differ.
        """
        if a != b:
            only_in_a = a - b
            only_in_b = b - a
            msg = f"{name} mismatch:\n"
            if only_in_a:
                msg += f"  Only in first: {only_in_a}\n"
            if only_in_b:
                msg += f"  Only in second: {only_in_b}"
            raise AssertionError(msg)
    
    @staticmethod
    def dataframe_equal(df1: pd.DataFrame, df2: pd.DataFrame,
                       check_dtype: bool = True, rtol: float = 1e-5):
        """
        Assert two DataFrames are equal, with detailed diff on failure.
        
        Parameters
        ----------
        df1, df2 : pd.DataFrame
            DataFrames to compare.
        check_dtype : bool
            Whether to check column dtypes.
        rtol : float
            Relative tolerance for numeric columns.
        
        Raises
        ------
        AssertionError
            If DataFrames differ.
        """
        try:
            pd.testing.assert_frame_equal(
                df1, df2, check_dtype=check_dtype, rtol=rtol, atol=1e-8
            )
        except AssertionError as e:
            raise AssertionError(str(e))
