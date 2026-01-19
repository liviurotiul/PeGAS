import sys

try:
    from .pipeline import parse_arguments, run_pipeline
except ImportError:
    import os
    package_root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    if package_root not in sys.path:
        sys.path.insert(0, package_root)
    from pegas.pipeline import parse_arguments, run_pipeline


def main(argv=None):
    args = parse_arguments(argv)
    run_pipeline(args)


if __name__ == "__main__":
    main(sys.argv[1:])
