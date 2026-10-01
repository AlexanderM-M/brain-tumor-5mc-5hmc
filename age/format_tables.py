"""Format publication-ready de-identified TSV tables; no model fitting."""
from pathlib import Path
import argparse
import csv
import re
from openpyxl import Workbook
from openpyxl.styles import Alignment,Font

COHORT=['Sample ID','Diagnosis','Age at sampling (years)','Sex','Sampling status',
        'DNAm age (years)','Apparent DNAm-age acceleration (years)','WGS-derived purity',
        'WGS purity status','Post-sampling follow-up (days)','Vital status','Survival-analysis inclusion']
SURVIVAL=['Section','Analysis','Population','Omitted sample ID','n','Deaths','Censored',
          'Median follow-up among censored patients (days)','Predictor','Scaling','HR',
          '95% CI lower','95% CI upper','Likelihood-ratio P','Direction relative to HR=1','HR range','Notes']


def read(path,columns):
    with path.open() as f:
        reader=csv.DictReader(f,delimiter='\t')
        if reader.fieldnames!=columns:raise ValueError('Unexpected columns; provide publication-ready tables only.')
        rows=list(reader)
    for row in rows:
        for value in row.values():
            if re.search(r'\bN\d{2}\.\d+\b|\b\d{4}-\d{2}-\d{2}\b|\b\d{1,2}[./]\d{1,2}[./]\d{4}\b',value):
                raise ValueError('Potential identifier/date in publication table.')
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort',type=Path,required=True);p.add_argument('--survival',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    tables=[read(a.cohort,COHORT),read(a.survival,SURVIVAL)]
    w=Workbook();w.remove(w.active);w.properties.creator='AGE manuscript authors';w.properties.lastModifiedBy='AGE manuscript preparation'
    numeric={'Age at sampling (years)','DNAm age (years)','Apparent DNAm-age acceleration (years)',
             'WGS-derived purity','Post-sampling follow-up (days)','n','Deaths','Censored',
             'Median follow-up among censored patients (days)','HR','95% CI lower','95% CI upper','Likelihood-ratio P'}
    for i,(columns,rows) in enumerate(zip([COHORT,SURVIVAL],tables),1):
        s=w.create_sheet('Supplementary Table '+str(i));s.append(columns)
        for r in rows:s.append([float(r[k]) if k in numeric and r[k]!='NA' else r[k] for k in columns])
        s.freeze_panes='A2';s.auto_filter.ref=s.dimensions;s.row_dimensions[1].height=58
        for row in s:
            for cell in row:
                cell.font=Font(name='Calibri',size=11,bold=cell.row==1)
                cell.alignment=Alignment(wrap_text=True,vertical='top')
                header=columns[cell.column-1]
                if header in ['HR','95% CI lower','95% CI upper','Likelihood-ratio P','WGS-derived purity']:
                    cell.number_format='0.000' if not isinstance(cell.value,(int,float)) or cell.value==0 or abs(cell.value)>=.001 else '0.000000'
                elif header in ['n','Deaths','Censored','Post-sampling follow-up (days)']:cell.number_format='0'
                elif header in numeric:cell.number_format='0.00' if 'years' in header else '0.#'
        for cell in s[1]:s.column_dimensions[cell.column_letter].width=min(42,max(15,len(cell.value)*.7))
        s.sheet_view.showGridLines=False;s.print_title_rows='1:1';s.page_setup.orientation='landscape'
        s.page_setup.paperSize=s.PAPERSIZE_A3;s.sheet_properties.pageSetUpPr.fitToPage=True;s.page_setup.fitToWidth=1;s.page_setup.fitToHeight=0
    a.output.parent.mkdir(parents=True,exist_ok=True);w.save(a.output)


if __name__=='__main__':main()
