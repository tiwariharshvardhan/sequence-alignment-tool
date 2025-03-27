import tkinter as tk
from tkinter import ttk, messagebox
from alignment_module import SequenceAligner

class GenomeAlignmentApp:
    def __init__(self, master):
        self.master = master
        master.title("Genome Sequence Alignment Tool")
        master.geometry("800x700")

        # Create aligner instance
        self.aligner = SequenceAligner()

        # Style configuration
        self.style = ttk.Style()
        self.style.configure("TLabel", font=("Arial", 10))
        self.style.configure("TButton", font=("Arial", 10))

        # Create main frame
        self.main_frame = ttk.Frame(master, padding="10")
        self.main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Sequence Input
        ttk.Label(self.main_frame, text="Sequence 1:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.seq1_entry = tk.Text(self.main_frame, height=5, width=80)
        self.seq1_entry.grid(row=1, column=0, columnspan=2, pady=5)

        ttk.Label(self.main_frame, text="Sequence 2:").grid(row=2, column=0, sticky=tk.W, pady=5)
        self.seq2_entry = tk.Text(self.main_frame, height=5, width=80)
        self.seq2_entry.grid(row=3, column=0, columnspan=2, pady=5)

        # Alignment Type Selection
        ttk.Label(self.main_frame, text="Alignment Type:").grid(row=4, column=0, sticky=tk.W, pady=5)
        self.alignment_type = tk.StringVar(value="global")
        self.global_radio = ttk.Radiobutton(
            self.main_frame, text="Global Alignment", 
            variable=self.alignment_type, value="global"
        )
        self.local_radio = ttk.Radiobutton(
            self.main_frame, text="Local Alignment", 
            variable=self.alignment_type, value="local"
        )
        self.global_radio.grid(row=5, column=0, sticky=tk.W)
        self.local_radio.grid(row=5, column=1, sticky=tk.W)

        # Alignment Parameters
        ttk.Label(self.main_frame, text="Match Score:").grid(row=6, column=0, sticky=tk.W, pady=5)
        self.match_score = tk.IntVar(value=2)
        self.match_entry = ttk.Entry(self.main_frame, textvariable=self.match_score, width=10)
        self.match_entry.grid(row=6, column=1, sticky=tk.W)

        ttk.Label(self.main_frame, text="Mismatch Penalty:").grid(row=7, column=0, sticky=tk.W, pady=5)
        self.mismatch_score = tk.IntVar(value=-1)
        self.mismatch_entry = ttk.Entry(self.main_frame, textvariable=self.mismatch_score, width=10)
        self.mismatch_entry.grid(row=7, column=1, sticky=tk.W)

        ttk.Label(self.main_frame, text="Gap Penalty:").grid(row=8, column=0, sticky=tk.W, pady=5)
        self.gap_penalty = tk.IntVar(value=-1)
        self.gap_entry = ttk.Entry(self.main_frame, textvariable=self.gap_penalty, width=10)
        self.gap_entry.grid(row=8, column=1, sticky=tk.W)

        # Align Button
        self.align_button = ttk.Button(self.main_frame, text="Align Sequences", command=self.perform_alignment)
        self.align_button.grid(row=9, column=0, columnspan=2, pady=10)

        # Results Display
        ttk.Label(self.main_frame, text="Alignment Results:").grid(row=10, column=0, sticky=tk.W, pady=5)
        self.results_text = tk.Text(self.main_frame, height=10, width=80, state='disabled')
        self.results_text.grid(row=11, column=0, columnspan=2, pady=5)

    def perform_alignment(self):
        # Validate input
        seq1 = self.seq1_entry.get("1.0", tk.END).strip().upper()
        seq2 = self.seq2_entry.get("1.0", tk.END).strip().upper()

        if not seq1 or not seq2:
            messagebox.showerror("Error", "Please enter both sequences")
            return

        # Update aligner parameters
        self.aligner = SequenceAligner(
            match_score=self.match_score.get(),
            mismatch_score=self.mismatch_score.get(),
            gap_penalty=self.gap_penalty.get()
        )

        # Perform alignment
        try:
            if self.alignment_type.get() == "global":
                result = self.aligner.needleman_wunsch(seq1, seq2)
                alignment_type = "Global Alignment"
            else:
                result = self.aligner.smith_waterman(seq1, seq2)
                alignment_type = "Local Alignment"

            # Display results
            self.results_text.config(state='normal')
            self.results_text.delete('1.0', tk.END)
            self.results_text.insert(tk.END, f"Alignment Type: {alignment_type}\n")
            self.results_text.insert(tk.END, f"Alignment Score: {result['score']}\n")
            self.results_text.insert(tk.END, f"Similarity: {result['similarity_percentage']:.2f}%\n\n")
            self.results_text.insert(tk.END, f"Aligned Sequence 1: {result['aligned_seq1']}\n")
            self.results_text.insert(tk.END, f"Aligned Sequence 2: {result['aligned_seq2']}")
            self.results_text.config(state='disabled')

        except Exception as e:
            messagebox.showerror("Alignment Error", str(e))

def main():
    root = tk.Tk()
    app = GenomeAlignmentApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()