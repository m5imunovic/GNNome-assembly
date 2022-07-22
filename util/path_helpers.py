from pathlib import Path


def get_project_root() -> Path:
    return Path(__file__).parent.parent


def get_vendor_path() -> Path:
    return get_project_root() / 'vendor'


def get_data_path() -> Path:
    return get_project_root() / 'data'


def get_ref_path() -> Path:
    return get_data_path() / 'references'