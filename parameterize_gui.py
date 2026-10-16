#!/usr/bin/env python3
# Author: Khatereh Kashmari
# Date: 2026-08-10
# Version: 9.1.0
# Required packages: Tkinter; dependencies of the selected parameterization backend
"""Desktop GUI for the unified LAMMPS parameterization dispatcher."""

from __future__ import annotations

import queue
import subprocess
import sys
import threading
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText


ROOT = Path(__file__).resolve().parent


class ParameterizationGUI(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("Atom Typing and Force Field Parameterization Tool for LAMMPS")
        self.geometry("940x760")
        self.minsize(800, 660)
        self.configure(background="#eef2f7")
        self.events: queue.Queue[tuple[str, object]] = queue.Queue()
        self.process: subprocess.Popen[str] | None = None

        self.input_var = tk.StringVar()
        self.ff_var = tk.StringVar(value="gaff2")
        self.output_var = tk.StringVar()
        self.charge_var = tk.StringVar(value="0")
        self.buffer_var = tk.StringVar(value="5.0")
        self.converter_var = tk.StringVar()
        self.frc_var = tk.StringVar()
        self.frc_charge_var = tk.StringVar(value="auto")
        self.use_existing_nta_var = tk.BooleanVar(value=False)
        self.nta_var = tk.StringVar()
        self.status_var = tk.StringVar(value="Ready")
        self._build()
        self._update_fields()
        self.after(100, self._drain_events)

    def _build(self) -> None:
        style = ttk.Style(self)
        style.theme_use("clam")
        style.configure("App.TFrame", background="#eef2f7")
        style.configure("Card.TLabelframe", background="#ffffff", borderwidth=1,
                        relief="solid", padding=14)
        style.configure("Card.TLabelframe.Label", background="#ffffff",
                        foreground="#19324d", font=("Helvetica", 12, "bold"))
        style.configure("Card.TLabel", background="#ffffff", foreground="#334e68")
        style.configure("Hint.TLabel", background="#ffffff", foreground="#627d98",
                        font=("Helvetica", 10))
        style.configure("AppTitle.TLabel", background="#eef2f7", foreground="#19324d",
                        font=("Helvetica", 13, "bold"))
        style.configure("AppStatus.TLabel", background="#eef2f7", foreground="#486581",
                        font=("Helvetica", 10, "bold"))
        style.configure("Footer.TLabel", background="#eef2f7", foreground="#829ab1",
                        font=("Helvetica", 9))
        style.configure("Accent.TButton", background="#1677ff", foreground="#ffffff",
                        font=("Helvetica", 11, "bold"), padding=(18, 9), borderwidth=0)
        style.map("Accent.TButton", background=[("active", "#0b63ce"),
                                                ("disabled", "#9fb3c8")])
        style.configure("Secondary.TButton", padding=(12, 7))
        style.configure("Status.Horizontal.TProgressbar", troughcolor="#d9e2ec",
                        background="#2cb67d", bordercolor="#d9e2ec", lightcolor="#2cb67d",
                        darkcolor="#2cb67d")

        header = tk.Frame(self, background="#102a43", height=104)
        header.pack(fill="x")
        header.pack_propagate(False)
        heading = tk.Frame(header, background="#102a43")
        heading.pack(fill="both", expand=True, padx=24, pady=17)
        tk.Label(heading, text="Atom Typing and Force Field Parameterization",
                 background="#102a43", foreground="#ffffff",
                 font=("Helvetica", 20, "bold")).pack(anchor="w")
        tk.Label(heading, text="LAMMPS data-file preparation and coefficient assignment",
                 background="#102a43", foreground="#bcccdc",
                 font=("Helvetica", 11)).pack(anchor="w", pady=(3, 0))

        frame = ttk.Frame(self, padding=(22, 18), style="App.TFrame")
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(1, weight=1)

        ttk.Label(frame, text="1  Select structure and force field", style="AppTitle.TLabel").grid(
            row=0, column=0, columnspan=3, sticky="w", pady=(0, 9))

        self._path_row(frame, 1, "Input structure", self.input_var, self._choose_input)

        ttk.Label(frame, text="Force field").grid(row=2, column=0, sticky="w", pady=7)
        ff_frame = ttk.Frame(frame)
        ff_frame.grid(row=2, column=1, columnspan=2, sticky="w")
        for text, value in (("GAFF2", "gaff2"), ("UFF", "uff"),
                            ("Existing .frc (auto type + build)", "frc")):
            ttk.Radiobutton(ff_frame, text=text, value=value, variable=self.ff_var,
                            command=self._update_fields).pack(side="left", padx=(0, 18))

        self._path_row(frame, 3, "Output directory", self.output_var, self._choose_output)

        self.gaff_frame = ttk.LabelFrame(frame, text="GAFF2 settings", style="Card.TLabelframe")
        self.gaff_frame.grid(row=4, column=0, columnspan=3, sticky="ew", pady=8)
        self.gaff_frame.columnconfigure(1, weight=1)
        ttk.Label(self.gaff_frame, text="Net charge").grid(row=0, column=0, sticky="w")
        ttk.Entry(self.gaff_frame, textvariable=self.charge_var, width=10).grid(
            row=0, column=1, sticky="w", padx=8)
        ttk.Label(self.gaff_frame, text="Box padding (Å)").grid(row=0, column=2, sticky="w")
        ttk.Entry(self.gaff_frame, textvariable=self.buffer_var, width=10).grid(
            row=0, column=3, sticky="w", padx=8)
        ttk.Label(self.gaff_frame, text="AMBER2LAMMPS (blank = auto)").grid(
            row=1, column=0, sticky="w", pady=(8, 0))
        ttk.Entry(self.gaff_frame, textvariable=self.converter_var).grid(
            row=1, column=1, columnspan=2, sticky="ew", padx=8, pady=(8, 0))
        ttk.Button(self.gaff_frame, text="Browse…", command=self._choose_converter).grid(
            row=1, column=3, pady=(8, 0))

        self.frc_frame = ttk.LabelFrame(frame, text=".frc parameterization settings",
                                        style="Card.TLabelframe")
        self.frc_frame.grid(row=5, column=0, columnspan=3, sticky="ew", pady=8)
        self.frc_frame.columnconfigure(1, weight=1)
        ttk.Label(self.frc_frame, text="Force-field .frc").grid(row=0, column=0, sticky="w")
        ttk.Entry(self.frc_frame, textvariable=self.frc_var).grid(row=0, column=1, sticky="ew", padx=8)
        ttk.Button(self.frc_frame, text="Browse…", command=self._choose_frc).grid(row=0, column=2)
        ttk.Label(self.frc_frame, text="Charge method").grid(row=1, column=0, sticky="w", pady=(8, 0))
        self.frc_charge_box = ttk.Combobox(
            self.frc_frame, textvariable=self.frc_charge_var, state="readonly",
            values=("auto", "gasteiger", "bond_increments", "fixed", "am1bcc", "input")
        )
        self.frc_charge_box.grid(row=1, column=1, sticky="w", padx=8, pady=(8, 0))
        ttk.Checkbutton(
            self.frc_frame, text="Use an existing .nta file (manual atom types)",
            variable=self.use_existing_nta_var, command=self._update_nta_mode
        ).grid(row=2, column=0, columnspan=3, sticky="w", pady=(10, 0))
        self.nta_frame = ttk.Frame(self.frc_frame)
        self.nta_frame.grid(row=3, column=0, columnspan=3, sticky="ew", pady=(8, 0))
        self.nta_frame.columnconfigure(1, weight=1)
        ttk.Label(self.nta_frame, text="Existing .nta").grid(row=0, column=0, sticky="w")
        ttk.Entry(self.nta_frame, textvariable=self.nta_var).grid(
            row=0, column=1, sticky="ew", padx=8)
        ttk.Button(self.nta_frame, text="Browse…", command=self._choose_nta).grid(row=0, column=2)
        self.nta_help = ttk.Label(self.frc_frame)
        self.nta_help.grid(row=4, column=0, columnspan=3, sticky="w", pady=(8, 0))

        button_frame = ttk.Frame(frame, style="App.TFrame")
        button_frame.grid(row=6, column=0, columnspan=3, sticky="ew", pady=(8, 8))
        self.run_button = ttk.Button(button_frame, text="▶  Run parameterization",
                                     command=self._start, style="Accent.TButton")
        self.run_button.pack(side="left")
        self.stop_button = ttk.Button(button_frame, text="Stop", command=self._stop,
                                      state="disabled", style="Secondary.TButton")
        self.stop_button.pack(side="left", padx=8)
        ttk.Label(button_frame, textvariable=self.status_var,
                  style="AppStatus.TLabel").pack(side="right")
        self.progress = ttk.Progressbar(button_frame, mode="indeterminate", length=160,
                                        style="Status.Horizontal.TProgressbar")
        self.progress.pack(side="right", padx=12)

        self.log = ScrolledText(frame, height=18, wrap="word", font=("Menlo", 10),
                                background="#0b1f33", foreground="#d9e2ec",
                                insertbackground="#ffffff", selectbackground="#276fbf",
                                relief="flat", padx=12, pady=10)
        self.log.grid(row=7, column=0, columnspan=3, sticky="nsew")
        frame.rowconfigure(7, weight=1)

        ttk.Label(frame, text="Atom Typing and Force Field Parameterization Tool for LAMMPS  •  Version 9.1.0  •  Khatereh Kashmari",
                  style="Footer.TLabel").grid(row=8, column=0, columnspan=3,
                                              sticky="e", pady=(8, 0))

    def _path_row(self, parent, row, label, variable, command) -> None:
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=7)
        ttk.Entry(parent, textvariable=variable).grid(row=row, column=1, sticky="ew", padx=8)
        ttk.Button(parent, text="Browse…", command=command).grid(row=row, column=2)

    def _choose_input(self) -> None:
        path = filedialog.askopenfilename(filetypes=[
            ("Supported structures", "*.mol2 *.vasp *.poscar *.data"), ("All files", "*")])
        if path:
            self.input_var.set(path)
            source = Path(path)
            self.output_var.set(str(source.parent / f"{source.stem}_{self.ff_var.get().upper()}"))

    def _choose_output(self) -> None:
        path = filedialog.askdirectory(mustexist=False)
        if path:
            self.output_var.set(path)

    def _choose_converter(self) -> None:
        path = filedialog.askopenfilename(filetypes=[("Python", "*.py"), ("All files", "*")])
        if path:
            self.converter_var.set(path)

    def _choose_frc(self) -> None:
        path = filedialog.askopenfilename(filetypes=[("Force field", "*.frc"), ("All files", "*")])
        if path:
            self.frc_var.set(path)

    def _choose_nta(self) -> None:
        path = filedialog.askopenfilename(filetypes=[("Atom types", "*.nta"), ("All files", "*")])
        if path:
            self.nta_var.set(path)

    def _update_nta_mode(self) -> None:
        manual = self.use_existing_nta_var.get()
        if manual:
            self.nta_frame.grid()
            self.frc_charge_box.configure(state="disabled")
            self.nta_help.configure(
                text="Manual mode: input must be a LAMMPS .data file; atom types come from the selected .nta.")
        else:
            self.nta_frame.grid_remove()
            self.frc_charge_box.configure(state="readonly")
            self.nta_help.configure(
                text="Automatic mode: atom typing creates the .nta before building the LAMMPS data file.")

    def _update_fields(self) -> None:
        if self.ff_var.get() == "gaff2":
            self.gaff_frame.grid()
        else:
            self.gaff_frame.grid_remove()
        if self.ff_var.get() == "frc":
            self.frc_frame.grid()
            self._update_nta_mode()
        else:
            self.frc_frame.grid_remove()
        if self.input_var.get():
            source = Path(self.input_var.get())
            self.output_var.set(str(source.parent / f"{source.stem}_{self.ff_var.get().upper()}"))

    def _command(self) -> list[str]:
        source = Path(self.input_var.get()).expanduser()
        if not source.is_file():
            raise ValueError("Select an existing input structure.")
        if not self.output_var.get().strip():
            raise ValueError("Select an output directory.")
        command = [sys.executable, str(ROOT / "parameterize.py"), str(source),
                   "-o", self.output_var.get().strip()]
        ff = self.ff_var.get()
        if ff == "gaff2":
            int(self.charge_var.get())
            float(self.buffer_var.get())
            command += ["--ff", "gaff2", "--charge", self.charge_var.get(),
                        "--buffer", self.buffer_var.get()]
            if self.converter_var.get().strip():
                command += ["--amber-converter", self.converter_var.get().strip()]
        elif ff == "uff":
            command += ["--ff", "uff"]
        else:
            frc = Path(self.frc_var.get())
            if not frc.is_file() or frc.suffix.lower() != ".frc":
                raise ValueError("Select an existing .frc file.")
            command += ["--ff", str(frc)]
            if self.use_existing_nta_var.get():
                nta = Path(self.nta_var.get()).expanduser()
                if not nta.is_file() or nta.suffix.lower() != ".nta":
                    raise ValueError("Select an existing .nta file.")
                if source.suffix.lower() != ".data":
                    raise ValueError("Manual .nta mode requires a LAMMPS .data input file.")
                command += ["--nta", str(nta)]
            else:
                command += ["--frc-charge", self.frc_charge_var.get()]
        return command

    def _start(self) -> None:
        try:
            command = self._command()
        except (ValueError, OSError) as exc:
            messagebox.showerror("Cannot start", str(exc))
            return
        self.log.delete("1.0", "end")
        self.log.insert("end", "$ " + " ".join(command) + "\n\n")
        self.status_var.set("Running…")
        self.progress.start(12)
        self.run_button.configure(state="disabled")
        self.stop_button.configure(state="normal")
        threading.Thread(target=self._worker, args=(command,), daemon=True).start()

    def _worker(self, command: list[str]) -> None:
        try:
            self.process = subprocess.Popen(command, cwd=ROOT, text=True,
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.STDOUT, bufsize=1)
            assert self.process.stdout is not None
            for line in self.process.stdout:
                self.events.put(("line", line))
            code = self.process.wait()
            self.events.put(("done", code))
        except Exception as exc:
            self.events.put(("error", str(exc)))

    def _stop(self) -> None:
        if self.process and self.process.poll() is None:
            self.process.terminate()
            self.status_var.set("Stopping…")

    def _drain_events(self) -> None:
        try:
            while True:
                kind, value = self.events.get_nowait()
                if kind == "line":
                    self.log.insert("end", str(value))
                    self.log.see("end")
                elif kind == "done":
                    code = int(value)
                    self.progress.stop()
                    if code == 0:
                        self.status_var.set("Completed")
                    elif code == 3:
                        self.status_var.set("Completed with missing coefficients")
                    else:
                        self.status_var.set(f"Failed (exit {code})")
                    self.run_button.configure(state="normal")
                    self.stop_button.configure(state="disabled")
                    if code == 0:
                        messagebox.showinfo("Complete", f"Files written to:\n{self.output_var.get()}")
                    elif code == 3:
                        messagebox.showwarning(
                            "Generated with warnings",
                            "The .nta and LAMMPS data files were generated, but some "
                            "coefficients are missing. Review audit.json before simulation.")
                elif kind == "error":
                    self.progress.stop()
                    self.status_var.set("Failed")
                    self.log.insert("end", f"\n[GUI ERROR] {value}\n")
                    self.run_button.configure(state="normal")
                    self.stop_button.configure(state="disabled")
        except queue.Empty:
            pass
        self.after(100, self._drain_events)


if __name__ == "__main__":
    ParameterizationGUI().mainloop()
