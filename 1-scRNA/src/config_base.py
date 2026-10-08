"""
Base configuration system for reproducible analysis pipelines.

Loads YAML config files, resolves paths, and provides a Config object that all
downstream modules can use. Configuration is the single source of truth; no
hardcoded values in analysis code.

Usage:
    cfg = load_config(path='configs/analysis.yaml')
    # Access parameters
    cfg['signature']['n_genes']  # dict-style
    cfg.path('processed')        # resolved path object
    cfg.stage('abundance')       # get stage-specific config
"""

import pathlib
import yaml
from typing import Union, Dict, Any


class Config:
    """Configuration loader with path resolution and stage access."""
    
    def __init__(self, config_dict: Dict[str, Any], config_path: Union[str, pathlib.Path]):
        """
        Initialize config.
        
        Parameters
        ----------
        config_dict : dict
            Parsed YAML configuration.
        config_path : str or Path
            Path to the config file (used for resolving relative paths).
        """
        self._config = config_dict
        self.config_path = pathlib.Path(config_path)
        self.root = self.config_path.parent.parent  # assume config/analysis.yaml
    
    def __getitem__(self, key: str) -> Any:
        """Dict-style access to config parameters."""
        return self._config[key]
    
    def get(self, key: str, default=None) -> Any:
        """Dict-style get with default."""
        return self._config.get(key, default)
    
    def path(self, name: str) -> pathlib.Path:
        """
        Resolve a named path relative to the project root.
        
        Parameters
        ----------
        name : str
            Key from config['paths'], e.g., 'processed', 'results', 'raw'.
        
        Returns
        -------
        pathlib.Path
            Absolute path, created if missing.
        """
        rel = self._config['paths'][name]
        p = self.root / rel
        p.mkdir(parents=True, exist_ok=True)
        return p
    
    def stage(self, stage_name: str) -> Dict[str, Any]:
        """
        Get stage-specific configuration.
        
        Parameters
        ----------
        stage_name : str
            Stage name, e.g., 'ingest', 'preprocess', 'signature'.
        
        Returns
        -------
        dict
            Configuration for that stage.
        """
        return self._config.get('stages', {}).get(stage_name, {})
    
    @property
    def seed(self) -> int:
        """Random seed for reproducibility."""
        return self._config.get('seed', 42)


def load_config(path: Union[str, pathlib.Path] = None) -> Config:
    """
    Load YAML configuration file.
    
    Parameters
    ----------
    path : str or Path, optional
        Path to config file. If None, searches for configs/analysis.yaml relative
        to the current working directory.
    
    Returns
    -------
    Config
        Configuration object.
    """
    if path is None:
        path = pathlib.Path.cwd() / 'configs' / 'analysis.yaml'
    
    path = pathlib.Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    
    with open(path) as f:
        config_dict = yaml.safe_load(f)
    
    return Config(config_dict, path)
