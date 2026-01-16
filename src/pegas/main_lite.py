import sys

try:
    from .pipeline import parse_arguments, run_pipeline
except ImportError:
    from pipeline import parse_arguments, run_pipeline


def main(argv=None):
    args = parse_arguments(argv)
    run_pipeline(args)


if __name__ == "__main__":
    main(sys.argv[1:])
