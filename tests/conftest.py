"""pytest 設定。RDKit が未インストールの場合は全テストをスキップ。"""
import pytest


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "rdkit: marks tests that require RDKit"
    )


def pytest_collection_modifyitems(config, items):
    try:
        import rdkit  # noqa: F401
    except ImportError:
        skip_rdkit = pytest.mark.skip(reason="RDKit not installed. Run: uv sync")
        for item in items:
            item.add_marker(skip_rdkit)
