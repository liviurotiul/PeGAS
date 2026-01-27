import json
import os
import sys


def launch_gui(config_filename="pegas_config.json"):
    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox
    except Exception as exc:
        print(f"[pegas] GUI unavailable: {exc}", file=sys.stderr)
        return None

    result = {"args": None}

    root = tk.Tk()
    root.title("PeGAS")
    root.resizable(False, False)

    default_gc = os.path.join(os.path.dirname(os.path.realpath(__file__)), "gc_content.json")

    data_var = tk.StringVar()
    output_var = tk.StringVar()
    cores_var = tk.StringVar()
    shovill_var = tk.StringVar()
    prokka_var = tk.StringVar()
    roary_var = tk.StringVar()
    gc_var = tk.StringVar(value=default_gc)
    overwrite_var = tk.BooleanVar(value=False)
    interactive_report_var = tk.BooleanVar(value=False)

    def apply_config(config_data):
        if not isinstance(config_data, dict):
            return
        data_dir = config_data.get("data")
        if data_dir:
            data_var.set(data_dir)
        cores = config_data.get("cores")
        if cores is not None:
            cores_var.set(str(cores))
        shovill = config_data.get("shovill_cpu_cores")
        if shovill is not None:
            shovill_var.set(str(shovill))
        prokka = config_data.get("prokka_cpu_cores")
        if prokka is not None:
            prokka_var.set(str(prokka))
        roary = config_data.get("roary_cpu_cores")
        if roary is not None:
            roary_var.set(str(roary))
        gc_path = config_data.get("gc")
        if gc_path:
            gc_var.set(gc_path)
        overwrite = config_data.get("overwrite")
        if overwrite is not None:
            overwrite_var.set(bool(overwrite))
        interactive_report = config_data.get("interactive_report")
        if interactive_report is None:
            interactive_report = config_data.get("html_report")
        if interactive_report is not None:
            interactive_report_var.set(bool(interactive_report))

    def load_config_for_output(output_dir):
        output_dir = output_dir.strip()
        if not output_dir:
            return
        output_dir = os.path.abspath(os.path.expanduser(output_dir))
        config_path = os.path.join(output_dir, config_filename)
        if not os.path.isfile(config_path):
            return
        try:
            with open(config_path, "r") as f:
                config_data = json.load(f)
        except Exception as exc:
            messagebox.showerror("Config error", f"Failed to load config: {exc}")
            return
        apply_config(config_data)

    def browse_data():
        folder = filedialog.askdirectory(title="Select input FASTQ directory")
        if folder:
            data_var.set(folder)

    def browse_output():
        folder = filedialog.askdirectory(title="Select output directory")
        if folder:
            output_var.set(folder)
            load_config_for_output(folder)

    def browse_gc():
        file_path = filedialog.askopenfilename(
            title="Select GC content JSON",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        if file_path:
            gc_var.set(file_path)

    def parse_int(value, label):
        value = value.strip()
        if not value:
            return None
        try:
            parsed = int(value)
        except ValueError:
            messagebox.showerror("Invalid value", f"{label} must be an integer.")
            return "invalid"
        if parsed <= 0:
            messagebox.showerror("Invalid value", f"{label} must be greater than zero.")
            return "invalid"
        return parsed

    def on_run():
        data_dir = data_var.get().strip()
        output_dir = output_var.get().strip()

        if not data_dir:
            messagebox.showerror("Missing value", "Input FASTQ directory is required.")
            return
        if not os.path.isdir(data_dir):
            messagebox.showerror("Invalid path", f"Input directory not found: {data_dir}")
            return
        if not output_dir:
            messagebox.showerror("Missing value", "Output directory is required.")
            return

        cores_val = parse_int(cores_var.get(), "Cores")
        if cores_val == "invalid":
            return
        shovill_val = parse_int(shovill_var.get(), "Shovill cores")
        if shovill_val == "invalid":
            return
        prokka_val = parse_int(prokka_var.get(), "Prokka cores")
        if prokka_val == "invalid":
            return
        roary_val = parse_int(roary_var.get(), "Roary cores")
        if roary_val == "invalid":
            return

        gc_path = gc_var.get().strip()
        if gc_path and not os.path.isfile(gc_path):
            messagebox.showerror("Invalid path", f"GC JSON file not found: {gc_path}")
            return

        args_list = ["-d", data_dir, "-o", output_dir]
        if cores_val is not None:
            args_list += ["-c", str(cores_val)]
        if shovill_val is not None:
            args_list += ["--shovill-cpu-cores", str(shovill_val)]
        if prokka_val is not None:
            args_list += ["--prokka-cpu-cores", str(prokka_val)]
        if roary_val is not None:
            args_list += ["--roary-cpu-cores", str(roary_val)]
        if gc_path:
            args_list += ["--gc", gc_path]
        if overwrite_var.get():
            args_list.append("--overwrite")
        if interactive_report_var.get():
            args_list.append("--interactive")

        result["args"] = args_list
        root.destroy()

    def on_cancel():
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_cancel)

    row = 0
    tk.Label(root, text="Input FASTQ directory").grid(row=row, column=0, sticky="w", padx=8, pady=4)
    tk.Entry(root, textvariable=data_var, width=50).grid(row=row, column=1, padx=8, pady=4)
    tk.Button(root, text="Browse", command=browse_data).grid(row=row, column=2, padx=8, pady=4)
    row += 1

    tk.Label(root, text="Output directory").grid(row=row, column=0, sticky="w", padx=8, pady=4)
    output_entry = tk.Entry(root, textvariable=output_var, width=50)
    output_entry.grid(row=row, column=1, padx=8, pady=4)
    output_entry.bind("<FocusOut>", lambda event: load_config_for_output(output_var.get()))
    tk.Button(root, text="Browse", command=browse_output).grid(row=row, column=2, padx=8, pady=4)
    row += 1

    tk.Label(root, text="Total cores (optional)").grid(row=row, column=0, sticky="w", padx=8, pady=4)
    tk.Entry(root, textvariable=cores_var, width=12).grid(row=row, column=1, sticky="w", padx=8, pady=4)
    row += 1

    tk.Label(root, text="Shovill cores (optional)").grid(row=row, column=0, sticky="w", padx=8, pady=4)
    tk.Entry(root, textvariable=shovill_var, width=12).grid(row=row, column=1, sticky="w", padx=8, pady=4)
    row += 1

    tk.Label(root, text="Prokka cores (optional)").grid(row=row, column=0, sticky="w", padx=8, pady=4)
    tk.Entry(root, textvariable=prokka_var, width=12).grid(row=row, column=1, sticky="w", padx=8, pady=4)
    row += 1

    tk.Label(root, text="Roary cores (optional)").grid(row=row, column=0, sticky="w", padx=8, pady=4)
    tk.Entry(root, textvariable=roary_var, width=12).grid(row=row, column=1, sticky="w", padx=8, pady=4)
    row += 1

    tk.Label(root, text="GC content JSON (optional)").grid(row=row, column=0, sticky="w", padx=8, pady=4)
    tk.Entry(root, textvariable=gc_var, width=50).grid(row=row, column=1, padx=8, pady=4)
    tk.Button(root, text="Browse", command=browse_gc).grid(row=row, column=2, padx=8, pady=4)
    row += 1

    tk.Checkbutton(root, text="Overwrite output directory", variable=overwrite_var).grid(row=row, column=1, sticky="w", padx=8, pady=6)
    row += 1

    tk.Checkbutton(root, text="Generate interactive HTML report", variable=interactive_report_var).grid(row=row, column=1, sticky="w", padx=8, pady=4)
    row += 1

    tk.Button(root, text="Run", command=on_run).grid(row=row, column=1, sticky="e", padx=8, pady=8)
    tk.Button(root, text="Cancel", command=on_cancel).grid(row=row, column=2, padx=8, pady=8)

    root.mainloop()
    return result["args"]
