# Genome Sequence Alignment Tool

## Project Structure
```
genome_alignment_project/
│
├── alignment_module.py      # Core alignment algorithms
├── alignment_gui.py         # Tkinter GUI application
└── requirements.txt         # Python dependencies
```

## Installation Steps

1. Create a virtual environment
```bash
python3 -m venv venv
source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
```

2. Install dependencies
```bash
pip install numpy
```

3. Running the Application
```bash
python alignment_gui.py
```

## Features
- Global Alignment (Needleman-Wunsch)
- Local Alignment (Smith-Waterman)
- Customizable alignment parameters
- Similarity percentage calculation
- User-friendly GUI

## Usage Instructions
1. Enter two DNA/protein sequences
2. Select alignment type (Global/Local)
3. Adjust scoring parameters if needed
4. Click "Align Sequences"
5. View alignment results in the output box

## Alignment Parameters
- Match Score: Reward for matching nucleotides
- Mismatch Penalty: Penalty for mismatched nucleotides
- Gap Penalty: Penalty for inserting gaps

## Example Sequences
- DNA: `ATCG`, `ATGG`
- Protein: `MVHLTPEEKSAVTALWGK`, `MVEWTDNAYRRK`

## Limitations
- Works best with short to medium-length sequences
- Computationally intensive for very long sequences