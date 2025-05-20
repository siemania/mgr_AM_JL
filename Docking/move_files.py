import shutil
import os

def move_file(source_path, destination_dir):
    destination_path = os.path.join(destination_dir, os.path.basename(source_path))
    shutil.move(source_path, destination_path)
    return destination_path

def move_dlg_xml_files(source_dir, destination_dir):
    for file in os.listdir(source_dir):
        if file.endswith('.dlg') or file.endswith('.xml'):
            move_file(os.path.join(source_dir, file), destination_dir)

if __name__ == '__main__':
    pass