"""Colour TRUE cells green and FALSE cells red in an Excel file.

Usage: python3 colour_true_false.py file.xlsx
"""
import sys
from openpyxl import load_workbook
from openpyxl.styles import PatternFill

GREEN = PatternFill("solid", fgColor="C6EFCE")
RED = PatternFill("solid", fgColor="FFC7CE")

path = sys.argv[1]
wb = load_workbook(path)
for ws in wb.worksheets:
    for row in ws.iter_rows():
        for cell in row:
            if cell.value in ("TRUE", True):
                cell.fill = GREEN
            elif cell.value in ("FALSE", False):
                cell.fill = RED
wb.save(path)
print(f"coloured {path}")
