import sys
from tqdm import tqdm

try:
    from .pipeline import parse_arguments, run_pipeline
    from .gui import launch_gui
    from .main import CONFIG_FILENAME
except ImportError:
    import os
    package_root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    if package_root not in sys.path:
        sys.path.insert(0, package_root)
    from pegas.pipeline import parse_arguments, run_pipeline
    from pegas.gui import launch_gui
    from pegas.main import CONFIG_FILENAME


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if len(argv) == 0:
        gui_args = launch_gui(CONFIG_FILENAME)
        if not gui_args:
            tqdm.write("[pegas] No arguments provided. Exiting.")
            return
        args = parse_arguments(gui_args, lite_mode=True)
    else:
        args = parse_arguments(argv, lite_mode=True)
    run_pipeline(args)


if __name__ == "__main__":
    main(sys.argv[1:])
