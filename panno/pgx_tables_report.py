#!/usr/bin/python
# -*- coding: UTF-8 -*-

import time, os, base64, csv
from itertools import *
from pathlib import Path

def save_to_csv(data, filename, output_dir):
    """Save data to .csv file"""

    os.makedirs(output_dir, exist_ok=True)
    
    filepath = os.path.join(output_dir, filename)
    
    if isinstance(data, dict):
        # For dicts (summary)
        with open(filepath, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            for category, drugs in data.items():
                writer.writerow([category] + drugs)
    elif hasattr(data, 'to_csv'):
        # For pandas DataFrame
        data.to_csv(filepath, encoding='utf-8', index=False)
    else:
        # For dictionary lists
        keys = data[0].keys() if data else []
        with open(filepath, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(data)

def csv_report(race, summary, prescribing_info, multi_var, single_var, phenotype_predict, clinical_anno, output_dir, sample_id):
    """Generate .csv"""
    
    # Creating a directory for the results
    output_dir = f"{os.path.dirname(output_dir)}/{os.path.basename(output_dir).split('.')[0]}"
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Summary
    save_to_csv(summary, f'{sample_id}_summary.csv', output_dir)
    
    # 2. Prescribing Info
    prescribing_info['PharmGKB Link'] = "https://www.pharmgkb.org/guidelineAnnotation/" + prescribing_info['PAID']
    save_to_csv(prescribing_info, f'{sample_id}_prescribing_info.csv', output_dir)
    
    # 3. Diplotype Detail - Multi-variant
    multi_var['Warning'] = ''
    multi_var.loc[multi_var['Gene'] == "CYP2B6", 'Warning'] = 'CYP2B6*29, CYP2B6*30 not considered'
    multi_var.loc[multi_var['Gene'] == "CYP2C19", 'Warning'] = 'CYP2C19*36, CYP2C19*37 not considered'
    multi_var.loc[multi_var['Gene'] == "CYP2D6", 'Warning'] = 'CYP2D6*5, CYP2D6*13, CYP2D6*61, CYP2D6*63, CYP2D6*68 and CNVs not considered'
    multi_var.loc[multi_var['Gene'] == "SLCO1B1", 'Warning'] = 'SLCO1B1*48, SLCO1B1*49 not considered'
    save_to_csv(multi_var, f'{sample_id}_diplotype_multi_variant.csv', output_dir)
    
    # 4. Diplotype Detail - Single-variant
    save_to_csv(single_var, f'{sample_id}_diplotype_single_variant.csv', output_dir)
    
    # 5. Phenotype Prediction
    # Сначала собираем все возможные категории фенотипов
    all_categories = set()
    for drug in phenotype_predict['Drug'].unique():
        drug_sub = phenotype_predict[phenotype_predict['Drug'] == drug]
        all_categories.update(drug_sub['PhenotypeCategory'].unique())
    
    # Создаем список всех возможных столбцов
    phenotype_fieldnames = ['Drug'] + sorted(all_categories)
    
    # Формируем данные
    phenotype_data = []
    for drug in phenotype_predict['Drug'].unique():
        drug_sub = phenotype_predict[phenotype_predict['Drug'] == drug]
        
        item = {'Drug': drug}
        for category in all_categories:
            category_data = drug_sub[drug_sub['PhenotypeCategory'] == category]
            if not category_data.empty:
                item[category] = category_data.iloc[0]['Prediction']
            else:
                item[category] = '-'
        
        phenotype_data.append(item)
    
    # Сохраняем с явным указанием всех полей
    with open(os.path.join(output_dir, f'{sample_id}_phenotype_prediction.csv'), 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=phenotype_fieldnames)
        writer.writeheader()
        writer.writerows(phenotype_data)
    
    # 6. Clinical Annotation
    clinical_data = []
    for _, row in clinical_anno.iterrows():
        clinical_data.append({
            'Drug': row['Drug'],
            'Category': row['PhenotypeCategory'],
            'Gene': row['Gene'],
            'Variant': row['Variant'],
            'Diplotype': row['Diplotype'],
            'Evidence Level': row['EvidenceLevel'],
            'Phenotype': row['PAnnoPhenotype'],
            'PharmGKB ID': row['CAID'],
            'PharmGKB Link': f"https://www.pharmgkb.org/clinicalAnnotation/{row['CAID']}"
        })
    save_to_csv(clinical_data, f'{sample_id}_clinical_annotation.csv', output_dir)

    # Report metadata
    metadata = {
        'Sample ID': sample_id,
        'Biogeographic Group': race,
        'Report Time': time.asctime(time.localtime(time.time())),
        'Tool': 'PAnno v0.3.1',
        'Reference': 'https://github.com/PreMedKB/PAnno',
        'Phenotype Warnings: Avoid use': ', '.join(summary['Avoid']) if summary['Avoid'] else 'N/A',
        'Phenotype Warnings: Not annotated': ', '.join(summary['NotInAnno']) if summary['NotInAnno'] else 'N/A'
    }
    save_to_csv([metadata], f'{sample_id}_metadata.csv', output_dir)