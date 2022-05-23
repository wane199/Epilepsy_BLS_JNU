import os

root = '/media/wane/Elements/EPTLE222/TLE'
command = '/media/wane/Elements/MRIneg/ROBEX/runROBEX.sh'

for subj in os.listdir(root):
    t1_path = os.path.join(root,subj,'T1')
    t1_file = os.path.join(t1_path,os.listdir(t1_path)[0])
    out_name = t1_file.split('.')[0]+"_brain.nii"
    os.system(command+" "+os.path.join(t1_path,t1_file)+" "+os.path.join(t1_path,out_name))