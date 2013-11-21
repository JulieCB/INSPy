
import numpy
import openpyxl

def read(filename):
    wb = openpyxl.load_workbook(filename)
    sheets = []
    for sheet in wb.worksheets:
        values = [[cell.value for cell in col] for col in sheet.columns]
        array = numpy.array(values, object).T
        sheets.append(array)
    return sheets
