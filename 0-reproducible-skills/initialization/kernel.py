import os
import pathlib

# Literal constant defining required directories
REQUIRED_DIRS_SPEC = (
    ("data", ("raw", "metadata", "processed")),
    ("src", ()),
    ("workflows", ()),
    ("configs", ()),
    ("results", ()),
    ("figures", ()),
    ("reports", ()),
    ("tests", ()),
    ("environment", ()),
)

def init_project_structure(path):
    """
    Initialize a project directory with the standard structure.
    
    Creates all required top-level directories and subdirectories.
    Idempotent: skips directories that already exist.
    
    Args:
        path: str or Path. Directory to initialize.
        
    Returns:
        None
        
    Raises:
        OSError: If directory creation fails.
    """
    path = pathlib.Path(path)
    
    try:
        # Create top-level directory
        path.mkdir(parents=True, exist_ok=True)
        print(f"✓ {path.name} initialized")
        
        # Create required directories
        for top_dir, subdirs in REQUIRED_DIRS_SPEC:
            top_path = path / top_dir
            top_path.mkdir(exist_ok=True)
            
            for subdir in subdirs:
                sub_path = top_path / subdir
                sub_path.mkdir(exist_ok=True)
        
        print(f"✓ Standard structure created")
        
    except Exception as e:
        raise OSError(f"Failed to initialize {path}: {e}")

def validate_project_structure(path):
    """
    Validate that a project has the correct directory structure.
    
    Args:
        path: str or Path. Directory to validate.
        
    Returns:
        list of str: Issues found (empty if valid).
    """
    path = pathlib.Path(path)
    issues = []
    
    # Check that root exists
    if not path.exists():
        return [f"Directory does not exist: {path}"]
    
    if not path.is_dir():
        return [f"Not a directory: {path}"]
    
    # Check required top-level directories
    found_dirs = set()
    for item in path.iterdir():
        if item.is_dir() and not item.name.startswith("."):
            found_dirs.add(item.name)
    
    required_top = {d[0] for d in REQUIRED_DIRS_SPEC}
    
    # Check for missing required directories
    missing = required_top - found_dirs
    for d in missing:
        issues.append(f"Missing directory: {d}/")
    
    # Check for unexpected directories
    unexpected = found_dirs - required_top
    for d in unexpected:
        issues.append(f"Unexpected directory: {d}/")
    
    # Check data/ subdirectories
    data_path = path / "data"
    if data_path.exists():
        data_subdirs = {d.name for d in data_path.iterdir() if d.is_dir() and not d.name.startswith(".")}
        required_data = {"raw", "metadata", "processed"}
        
        missing_data = required_data - data_subdirs
        for d in missing_data:
            issues.append(f"data/ missing subdirectory: {d}/")
        
        unexpected_data = data_subdirs - required_data
        for d in unexpected_data:
            issues.append(f"data/ unexpected subdirectory: {d}/")
    
    return issues
