import pandas as pd
import numpy as np

DATA_DIR = "../../data/"
COVID_METADATA_XLSX = DATA_DIR + "/meta/covid_metadata.xlsx"


def load_meta_data():
    dfx = pd.read_excel(COVID_METADATA_XLSX,
                        index_col="sampleID")  # type: ignore
    filter_types = ["fresh BALF", "fresh Sputum", "fresh PFMC"]
    dfx = dfx.loc[~dfx["Sample type"].isin(filter_types),]
    vcx = dfx["Sample type"].value_counts()
    qcx = pd.read_excel(COVID_METADATA_XLSX, sheet_name="qc",
                        index_col="sampleID")  # type:ignore
    qcx = qcx.loc[dfx.index]
    qcx = qcx.loc[qcx["Ncells_manual_removal"] > 1000]
    dfx = dfx.loc[qcx.index]
    return dfx, qcx, vcx


def threshold_stats(qcx):
    for x in [200, 500, 800, 1000, 2000, 5000]:
        qcx3 = qcx.loc[qcx["Ncells_manual_removal"] > x]
        print("%4d %4d %8d" % (x, qcx3.shape[0],
                               sum(qcx3["Ncells_manual_removal"])))


def select_diseased_from(dfx, qcx, ndis, from_set):
    # disease_patients = list(dfx[dfx["SARS-CoV-2"] == "positive"].index)
    disease_patients = from_set
    disx = np.random.choice(range(len(disease_patients)), ndis, replace=False)
    disx = [disease_patients[x] for x in disx]
    nsum = sum(qcx.loc[disx].Ncells_manual_removal)
    print(nsum)
    return disx, nsum


def select_diseased(dfx, qcx, ndis):
    disease_patients = list(dfx[dfx["SARS-CoV-2"] == "positive"].index)
    disx = np.random.choice(range(len(disease_patients)), ndis, replace=False)
    disx = [disease_patients[x] for x in disx]
    nsum = sum(qcx.loc[disx].Ncells_manual_removal)
    print(nsum)
    return disx, nsum


def select_control_full(dfx, qcx):
    control_patients = list(dfx[dfx["SARS-CoV-2"] == "negative"].index)
    ctrlx = range(len(control_patients))
    ctrlx = [control_patients[x] for x in ctrlx]
    nsum = sum(qcx.loc[ctrlx].Ncells_manual_removal)
    print(nsum)
    return ctrlx, nsum


def select_control(dfx, qcx, ndis):
    control_patients = list(dfx[dfx["SARS-CoV-2"] == "negative"].index)
    ctrlx = np.random.choice(range(len(control_patients)), ndis, replace=False)
    ctrlx = [control_patients[x] for x in ctrlx]
    nsum = sum(qcx.loc[ctrlx].Ncells_manual_removal)
    print(nsum)
    return ctrlx, nsum


def ctrl50kp():
    return ["S-HC025", "S-HC010", "S-HC013", "S-HC008", "S-HC019-2"]


def dis10kp():
    return ['S-S062', 'S-M010-3', 'S-S054']


def dis25Kp():
    return ['S-S035-1', 'S-S037', 'S-S013', 'S-M077', 'S-S036-3']


def dis150kp():
    return ['S-M061-2', 'S-M041-1', 'S-S021-4', 'S-M044-1', 'S-M043-2',
            'S-M053', 'S-S018', 'S-S036-3', 'S-S052', 'S-S035-3',
            'S-M056', 'S-S087-2', 'S-M049', 'S-M020', 'S-M001',
            'S-S082', 'S-M035-1', 'S-M012', 'S-S083', 'S-S050',
            'S-S027', 'S-M018', 'S-S086-2', 'S-S061']


def dis50kp():
    return ["S-M035-1", "S-S035-2", "S-M061-2", "S-M028", "S-S079"]


def dis500kp():
    return [
        "S-S068", "S-M018", "S-S050", "S-M066", "S-S013", "S-S035-1", "S-S001-2",
        "S-S081", "S-M037", "S-M042-2", "S-S078", "S-S047", "S-M048", "S-S018",
        "S-M042-1", "S-S050", "S-M010-2", "S-M009-2", "S-S084", "S-S044", "S-S047",
        "S-M073", "S-M050", "S-M016", "S-M008-2", "S-M078", "S-M011", "S-M054",
        "S-S092", "S-M007-4", "S-M064", "S-S082", "S-M004-1", "S-S068", "S-S035-4",
        "S-M030", "S-S013", "S-M004-6", "S-S092", "S-M033", "S-S057", "S-M041-1",
        "S-M063-1", "S-M010-1", "S-M069", "S-M049", "S-M059-2", "S-S017", "S-S082",
        "S-M005",
    ]


def ds50():
    return [
        'S-M018', 'S-S050', 'S-M066', 'S-S013', 'S-S035-1',
        'S-S001-2', 'S-S081', 'S-M037'
    ]


def ds350_r1():
    return [
        'S-M045', 'S-S074-1', 'S-M042-2', 'S-S022-1', 'S-S063', 'S-M044-1',
        'S-S036-3', 'S-S036-3', 'S-M074-2', 'S-S087-2', 'S-M036-1', 'S-M067',
        'S-M024', 'S-M007-2', 'S-S036-3', 'S-M077', 'S-M010-6', 'S-M004-2',
        'S-M034', 'S-M055', 'S-S039', 'S-S032-3', 'S-M059-2', 'S-S073-1', 'S-S081',
        'S-M070', 'S-S091', 'S-S022-5', 'S-M004-2', 'S-S091', 'S-M028', 'S-S050',
        'S-S022-2', 'S-M061-2', 'S-M018', 'S-S037', 'S-M036-2', 'S-S050', 'S-M038',
        'S-M046', 'S-M028', 'S-M034', 'S-S086-2', 'S-M007-5', 'S-M014', 'S-M021',
        'S-M063-1', 'S-S081', 'S-M004-4', 'S-M008-1', 'S-S021-5', 'S-M016',
        'S-M079', 'S-M009-3', 'S-S014'
    ]


def ds350():
    return [
        "S-S041", "S-S049", "S-M060-1", "S-S076-1", "S-M011", "S-M010-4", "S-S080",
        "S-M051", "S-S020", "S-S013", "S-S022-2", "S-S039", "S-M018", "S-M007-2",
        "S-M027", "S-M004-6", "S-M033", "S-M014", "S-S018", "S-S026", "S-S086-2",
        "S-S031", "S-M042-2", "S-S073-1", "S-M008-2", "S-S083", "S-S021-4", "S-S043",
        "S-M010-3", "S-S077", "S-M004-3", "S-M017", "S-S021-2", "S-M005", "S-M004-2",
        "S-M058-1", "S-S036-1", "S-M056", "S-S091", "S-S070-2", "S-M007-4",
        "S-M010-2", "S-M076-2", "S-M043-1", "S-M028", "S-S030", "S-S001-2",
        "S-S023", "S-S035-4", "S-M041-2", "S-M007-6", "S-S021-3", "S-S085-2",
        "S-S046", "S-M008-1", "S-S033", "S-M040-2", "S-M077", "S-S056", "S-M009-6"
    ]


def ds500_r1():
    return [
        'S-S077', 'S-M010-3', 'S-M011', 'S-S045', 'S-M009-5', 'S-S045',
        'S-M019', 'S-M066', 'S-S077', 'S-M039-2', 'S-M031-1', 'S-M036-1',
        'S-S073-1', 'S-M029', 'S-M062-2', 'S-M032', 'S-M008-2', 'S-M059-1',
        'S-S022-5', 'S-M016', 'S-S014', 'S-M022', 'S-M060-2', 'S-M048',
        'S-S088-2', 'S-M066', 'S-M071', 'S-S074-2', 'S-S017', 'S-M051',
        'S-M040-1', 'S-M032', 'S-M010-6', 'S-M072', 'S-M004-4', 'S-S070-3',
        'S-S020', 'S-S091', 'S-M035-1', 'S-M076-2', 'S-S084', 'S-S076-1',
        'S-M010-1', 'S-S035-2', 'S-S025', 'S-M062-2', 'S-S022-1', 'S-M006',
        'S-S045', 'S-M054', 'S-S067', 'S-M009-5', 'S-S041', 'S-S029',
        'S-M009-2', 'S-S066', 'S-S051', 'S-S057', 'S-M048', 'S-M061-2',
        'S-M019', 'S-S035-1', 'S-S039', 'S-M040-1', 'S-M037', 'S-S021-3',
        'S-M018', 'S-M058-1', 'S-S090-2', 'S-S061', 'S-M001', 'S-M071',
        'S-S064', 'S-M040-1', 'S-M054', 'S-M005', 'S-S089-2', 'S-S086-2',
        'S-M056', 'S-M026-1', 'S-S078', 'S-S034', 'S-M032', 'S-M060-2',
        'S-S087-2'
    ]


def ds500():
    return [
        "S-M043-1", "S-M071", "S-M004-3", "S-M079", "S-M007-1", "S-M055", "S-S087-2",
        "S-M014", "S-M077", "S-M051", "S-S015", "S-S035-4", "S-S047", "S-S085-2",
        "S-S049", "S-M026-1", "S-S040", "S-S083", "S-M072", "S-M053", "S-M035-2",
        "S-M026-3", "S-M009-4", "S-M016", "S-S034", "S-S063", "S-S051", "S-S031",
        "S-M048", "S-S025", "S-M004-6", "S-S026", "S-M043-2", "S-S054", "S-S036-3",
        "S-M015", "S-S084", "S-M066", "S-S090-2", "S-S020", "S-S080", "S-S022-3",
        "S-S021-4", "S-M059-1", "S-M007-3", "S-M007-5", "S-S077", "S-M004-2",
        "S-M009-5", "S-M017", "S-M063-1", "S-M044-1", "S-M037", "S-S023",
        "S-S070-3", "S-S076-2", "S-S013", "S-M010-4", "S-M059-2", "S-S024",
        "S-M041-1", "S-M039-1", "S-S027", "S-S043", "S-M040-1", "S-M026-2",
        "S-S037", "S-M067", "S-S028", "S-M004-5", "S-M042-2", "S-M036-1",
        "S-M062-2", "S-M004-1", "S-M021", "S-M022", "S-M042-1", "S-M039-2",
        "S-M018", "S-M078"
    ]


def ds700_r1():
    return [
        'S-M013', 'S-S018', 'S-M059-1', 'S-M007-5', 'S-M059-1', 'S-S041',
        'S-M004-1', 'S-M015', 'S-S029', 'S-S090-2', 'S-M076-2', 'S-S022-5',
        'S-S036-1', 'S-S091', 'S-M001', 'S-M007-5', 'S-S092', 'S-M010-6',
        'S-S064', 'S-M043-2', 'S-S018', 'S-S076-2', 'S-M044-2', 'S-S022-1',
        'S-S022-1', 'S-S021-2', 'S-S038', 'S-M010-6', 'S-M023', 'S-M070',
        'S-M026-3', 'S-M007-3', 'S-M032', 'S-M028', 'S-S068', 'S-M020',
        'S-S067', 'S-M064', 'S-S036-3', 'S-M013', 'S-S053', 'S-M058-1',
        'S-M026-3', 'S-M027', 'S-S081', 'S-S023', 'S-S022-2', 'S-S031', 'S-S030',
        'S-M056', 'S-M025', 'S-S026', 'S-M068', 'S-S062', 'S-M016', 'S-S052',
        'S-M015', 'S-M009-2', 'S-S022-3', 'S-S014', 'S-M010-4', 'S-M038',
        'S-M004-1', 'S-M074-2', 'S-M010-3', 'S-M008-1', 'S-M040-2',
        'S-M009-6', 'S-S047', 'S-S015', 'S-M077', 'S-S087-2', 'S-S018',
        'S-S016', 'S-M043-2', 'S-S014', 'S-S070-3', 'S-M007-3', 'S-S014',
        'S-M044-1', 'S-S081', 'S-S053', 'S-S064', 'S-S045', 'S-M054',
        'S-S027', 'S-M007-6', 'S-M063-2', 'S-S021-4', 'S-S013', 'S-S062',
        'S-M068', 'S-S063', 'S-M007-4', 'S-M048', 'S-M047', 'S-S089-2',
        'S-S021-1', 'S-M079', 'S-S032-3', 'S-M032', 'S-M040-2', 'S-S017',
        'S-S026', 'S-M030', 'S-S092', 'S-M012', 'S-M061-2', 'S-S036-3',
        'S-M028', 'S-M040-1', 'S-M007-5', 'S-M076-2', 'S-M068', 'S-M035-2',
        'S-M017', 'S-S075-1', 'S-M016', 'S-S035-2', 'S-M036-2'
    ]


def ds700():
    return [
        "S-S022-2", "S-S044", "S-S015", "S-M020", "S-M056", "S-M036-1",
        "S-S037", "S-S013", "S-S088-2", "S-M048", "S-M007-2", "S-M068",
        "S-M025", "S-S021-5", "S-S080", "S-S090-2", "S-M007-1", "S-S032-3",
        "S-S079", "S-S046", "S-S022-4", "S-S036-2", "S-M011", "S-M043-2",
        "S-M009-6", "S-M077", "S-M004-4", "S-S050", "S-M037", "S-S089-2",
        "S-M008-1", "S-S078", "S-M023", "S-M029", "S-S043", "S-M061-2",
        "S-S057", "S-M060-1", "S-M061-1", "S-M009-3", "S-S087-2", "S-M035-2",
        "S-S066", "S-S045", "S-M010-6", "S-S091", "S-S035-2", "S-M041-1",
        "S-S021-2", "S-M063-1", "S-M021", "S-S022-3", "S-S086-2", "S-M047",
        "S-M074-2", "S-S040", "S-S027", "S-M034", "S-M067", "S-S063",
        "S-S065", "S-S036-3", "S-M018", "S-S029", "S-M001", "S-S023",
        "S-S068", "S-M007-6", "S-M022", "S-S022-5", "S-S075-1", "S-M009-4",
        "S-M059-1", "S-S069-3", "S-S026", "S-S051", "S-M007-4", "S-M015",
        "S-S021-4", "S-M004-1", "S-M010-5", "S-M039-2", "S-M013", "S-S016",
        "S-M071", "S-S048", "S-S083", "S-M009-5", "S-M031-2", "S-M042-2",
        "S-M033", "S-S014", "S-M045", "S-M066", "S-M049", "S-M007-3",
        "S-S025", "S-M009-1", "S-M040-2", "S-S067", "S-S042", "S-M054",
        "S-M006", "S-M024", "S-M007-5", "S-S054", "S-S039", "S-M030",
        "S-M038", "S-M079", "S-S074-1", "S-S028", "S-S031", "S-M064", "S-S084"
    ]


def dsfull():
    return [
        'S-S070-2', 'S-S070-3', 'S-S069-3', 'S-M056', 'S-M044-1', 'S-M043-1',
        'S-M048', 'S-M044-2', 'S-M043-2', 'S-S054', 'S-S056', 'S-M042-1', 'S-M041-1',
        'S-M049', 'S-M046', 'S-M047', 'S-S055', 'S-S057', 'S-M045', 'S-M041-2',
        'S-M042-2', 'S-M055', 'S-M053', 'S-S067', 'S-S065', 'S-M051', 'S-S064',
        'S-M054', 'S-M052', 'S-S068', 'S-S066', 'S-S059', 'S-S060', 'S-M050', 'S-S061',
        'S-S062', 'S-S063', 'S-M061-1', 'S-M061-2', 'S-S073-1', 'S-S074-1', 'S-S074-2',
        'S-M062-2', 'S-M058-1', 'S-M058-2', 'S-M063-1', 'S-M063-2', 'S-M059-1',
        'S-M059-2', 'S-M060-1', 'S-M060-2', 'S-S075-1', 'S-S076-1', 'S-S076-2',
        'S-S035-1', 'S-S035-2', 'S-S035-3', 'S-S035-4', 'S-S036-1', 'S-S036-2',
        'S-S036-3', 'S-M064', 'S-M066', 'S-M067', 'S-M068', 'S-S091', 'S-S090-2',
        'S-S092', 'S-M076-2', 'S-M074-2', 'S-S088-2', 'S-S085-2', 'S-S089-2',
        'S-S087-2', 'S-S086-2', 'S-M077', 'S-M078', 'S-M079', 'S-S024', 'S-S026',
        'S-M013', 'S-M014', 'S-M015', 'S-M016', 'S-S023', 'S-S025', 'S-S027',
        'S-S034', 'S-M025', 'S-S033', 'S-M028', 'S-M029', 'S-M027', 'S-M026-1',
        'S-S032-3', 'S-M026-2', 'S-M026-3', 'S-M004-1', 'S-M004-2', 'S-M004-3',
        'S-S022-1', 'S-S022-2', 'S-S022-3', 'S-S022-4', 'S-S022-5', 'S-M010-1',
        'S-M010-2', 'S-M010-3', 'S-M010-4', 'S-M010-5', 'S-M010-6', 'S-M009-1',
        'S-M009-2', 'S-M011', 'S-M012', 'S-M005', 'S-M006', 'S-M004-4', 'S-M004-5',
        'S-M004-6', 'S-S013', 'S-S014', 'S-S015', 'S-S016', 'S-S017', 'S-S018',
        'S-S019', 'S-S020', 'S-M007-1', 'S-M007-2', 'S-M007-3', 'S-M007-4', 'S-M007-5',
        'S-M007-6', 'S-M008-1', 'S-M008-2', 'S-M009-3', 'S-M009-4', 'S-M009-5',
        'S-M009-6', 'S-S021-1', 'S-S021-2', 'S-S021-3', 'S-S021-4', 'S-S021-5',
        'S-M018', 'S-S029', 'S-S030', 'S-M023', 'S-S031', 'S-M024', 'S-M017', 'S-M019',
        'S-M020', 'S-M021', 'S-S028', 'S-M022', 'S-M001', 'S-M030', 'S-M031-1',
        'S-M031-2', 'S-M032', 'S-M033', 'S-M034', 'S-M035-1', 'S-M035-2', 'S-M036-1',
        'S-M036-2', 'S-M037', 'S-M038', 'S-M039-1', 'S-M039-2', 'S-M040-1', 'S-M040-2',
        'S-S001-2', 'S-M069', 'S-M070', 'S-M071', 'S-M072', 'S-M073', 'S-S077',
        'S-S078', 'S-S079', 'S-S080', 'S-S081', 'S-S082', 'S-S083', 'S-S084', 'S-S037',
        'S-S038', 'S-S039', 'S-S040', 'S-S041', 'S-S042', 'S-S043', 'S-S044', 'S-S045',
        'S-S046', 'S-S047', 'S-S048', 'S-S049', 'S-S050', 'S-S051', 'S-S052', 'S-S053'
    ]


def subset150_25k():
    return ['S-M053', 'S-S018', 'S-S086-2']


def subset150_50k():
    return [
        'S-M012',
        'S-M053',
        'S-M061-2',
        'S-S018',
        'S-S086-2',
        'S-S087-2',
    ]


def subset150_60k():
    return [
        'S-M012',
        'S-M044-1',
        'S-M053',
        'S-M061-2',
        'S-S018',
        'S-S082',
        'S-S086-2',
        'S-S087-2',
    ]


def subset150_75k():
    return [
        'S-M012',
        'S-M035-1',
        'S-M043-2',
        'S-M044-1',
        'S-M053',
        'S-M061-2',
        'S-S018',
        'S-S082',
        'S-S086-2',
        'S-S087-2',
    ]


def subset150_80k():
    return [
        'S-M012',
        'S-M035-1',
        'S-M043-2',
        'S-M044-1',
        'S-M049',
        'S-M053',
        'S-M061-2',
        'S-S018',
        'S-S082',
        'S-S086-2',
        'S-S087-2',
    ]


def subset150_100k():
    return [
        'S-M001',
        'S-M012',
        'S-M035-1',
        'S-M043-2',
        'S-M044-1',
        'S-M049',
        'S-M053',
        'S-M056',
        'S-M061-2',
        'S-S018',
        'S-S050',
        'S-S082',
        'S-S086-2',
        'S-S087-2',
    ]


def subset150_125k():
    return [
        'S-M001',
        'S-M012',
        'S-M020',
        'S-M035-1',
        'S-M043-2',
        'S-M044-1',
        'S-M049',
        'S-M053',
        'S-M056',
        'S-M061-2',
        'S-S018',
        'S-S027',
        'S-S035-3',
        'S-S036-3',
        'S-S050',
        'S-S082',
        'S-S083',
        'S-S086-2',
        'S-S087-2',
    ]


def get_dsfull():
    dfx, _, _ = load_meta_data()
    return list(dfx[dfx['SARS-CoV-2'] == 'positive'].index)   # type:ignore


def get_full():
    dfx, _, _ = load_meta_data()
    return list(dfx.index)   # type:ignore


def select_random_control(nsamples):
    dfx, qcx, _ = load_meta_data()
    ctrl_list, ctrl_sum = select_control(dfx, qcx, nsamples)
    return ctrl_list, ctrl_sum


def select_random_diseased(nsamples):
    dfx, qcx, _ = load_meta_data()
    diseased_list, diseased_sum = select_diseased(dfx, qcx, nsamples)
    return diseased_list, diseased_sum
