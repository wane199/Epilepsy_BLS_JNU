#coding=utf-8
# By rengaosheng
# email: iamrgs@foxmail.com
# update 2021年2月28日, 新增 modality 类别筛选分类

import os
import shutil
import pydicom


def is_dicom_file(filename):

    # 判断某文件是否是dicom格式的文件
    file_stream = open(filename, 'rb')
    file_stream.seek(128)
    data = file_stream.read(4)
    file_stream.close()
    if data == b'DICM':
        return True
    return False


def parseDicom(file_path, show_log=False, file_type=["MR", "CT", "PT", "SR"]):
    ds = pydicom.dcmread(file_path)
    modality = ds.Modality
    if modality not in file_type:
        print(f"{modality} is not in selected type: {file_type}. Pass.")
        return None
    first_dir = f"{ds.PatientName}_{ds.PatientAge}_{ds.PatientSex}" 
    second_dir = f"{ds.SeriesDescription}"
    #first_dir = f"HY{ds.StudyDate}_monkey_head_{ds.Modality}"
    #first_dir = f"HY{ds.StudyDate}_monkey_head_MRI"
    #second_dir = f"{ds.PatientName}_{ds.PatientAge}_{ds.PatientSex}"
    #second_dir = f"{ds.PatientName}_{ds.PatientSex}"
    #third_path = f"{ds.SeriesDescription}"

    #final_dir = f"{first_dir}/{second_dir}/{third_path}"
    final_dir = f"{first_dir}/{second_dir}"
    final_dir = final_dir.replace(" ", "_")       # replace the special symbol.
    final_dir = final_dir.replace("*", "_")
    final_dir = final_dir.replace("-", "_")
    final_dir = final_dir.replace(":", "_")

    if show_log:
        print("--------------------")
        print(f"File: {file_path}")
        print(f"StudyDate: {ds.StudyDate}")
        print(f"Modality: {modality}")

        print(f"PatientName: {ds.PatientName}")
        #print(f"PatientAge: {ds.PatientAge}")
        print(f"PatientSex: {ds.PatientSex}")

        print(f"SeriesDescription: {ds.SeriesDescription}")
        print("--------------------")

    return final_dir


def sortDicom(src_dicom_path, det_dicom_path, judge_dicom=True, show_log=False, file_type=["MR", "CT", "PT", "SR"]):
    '''
    src_dicom_path: The dicom dir to be sorted. EG: "C:/Users/wane/Desktop/monkey/Lijing/0807/"
    det_dicom_path: The dir to store the sorted dicom files.
    judge_dicom: judge the file is dicom or not. Default is True.
    file_type: The Modality type you want to process.
    '''
    # 1. search dicom files.
    print(f"Searching for dicom files...")
    det_dicom_path = det_dicom_path.replace("\\", "/")
    listFilesDCM  = []
    for dirName, subdirList, fileList in os.walk(src_dicom_path):
        for filename in fileList:
            file_path = os.path.join(dirName, filename)
            file_path = file_path.replace("\\", "/")
            if judge_dicom:
                is_dicom = is_dicom_file(file_path)
                if is_dicom is False:
                    print(f"{file_path} is not a dicom file. pass.")
                    continue
            listFilesDCM.append(file_path)
    total_num = len(listFilesDCM)
    print(f"Total: {total_num} dicom files. ")

    # 2. parse dicom file, copy into dir.
    file_idx = 1
    for file_path in listFilesDCM:
        final_dir = parseDicom(file_path, show_log=show_log, file_type=file_type)
        if final_dir is None:
            print(f"Current: {file_idx} / {total_num}. pass")
            file_idx += 1
            continue
        det_path = os.path.join(det_dicom_path, final_dir)
        if not os.path.exists(det_path):
            os.makedirs(det_path)
        shutil.copy(file_path, det_path)
        print(f"Current: {file_idx} / {total_num}")
        file_idx += 1
    print(f"------- Done. ------")


if __name__ == '__main__':

    src = "/media/wane/wade/ZHL/1/图片数据/2014年（69个）/2014-04-29 吴丽娇/PA0"
    det = "/media/wane/wade/ZHL/1/图片数据/2014年（69个）/2014-04-29 吴丽娇/PA0"
    sortDicom(src, det, show_log=True, file_type=["PT","MR","CT"])
    # sortDicom(src, det)