from pathlib import Path

ROOT_DIR = Path(__file__).parent.parent.parent
DATA_DIR = ROOT_DIR / 'data'
DEFAULT_PARAMETERS_DIR = ROOT_DIR / 'default_parameters'
PROJECT_DIR = ROOT_DIR / 'projects'
OUTPUT_DIR = ROOT_DIR / 'output'
DEFAULT_EXCLUDELISTS_DIR = DEFAULT_PARAMETERS_DIR / 'excludelists'

if __name__ == '__main__':
    print(__file__)
    print(ROOT_DIR)
    print(DATA_DIR)
    print(DEFAULT_PARAMETERS_DIR)
    print(PROJECT_DIR)