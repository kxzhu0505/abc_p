import os
import re
import sys


#将划分后的网表merged为一个新的网表，使用方法：./IncrMerg.py /需要划分的网表所在路径 /划分后网表保存路径

def process_blif_files(directory):
    for filename in os.listdir(directory):
        filepath = os.path.join(directory, filename)

        if os.path.isdir(filepath):
            process_blif_files(filepath)
        elif filename.endswith(".blif"):
            process_blif_file(filepath)

#将abc输出的blif文件的输入输出格式进行整理
def process_blif_file(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()

    inputs = set()  
    is_input_section = False
    for line in lines:
        line = line.strip()
        if line.startswith('.inputs'):
            is_input_section = True
            line = line.split('.inputs')[1].strip()
        if is_input_section:
            if line.endswith('\\'):
                line = line[:-1].strip()
            else:
                is_input_section = False
            inputs.update(line.split()) 

    outputs = set()
    is_output_section = False
    for line in lines:
        line = line.strip()
        if line.startswith('.outputs'):
            is_output_section = True
            line = line.split('.outputs')[1].strip()
        if is_output_section:
            if line.endswith('\\'):
                line = line[:-1].strip()
            else:
                is_output_section = False
            outputs.update(line.split())


    start_index = None
    for i, line in enumerate(lines):
        if line.startswith('.names'):
            start_index = i
            break


    new_lines = []
    new_lines.append('.model netlist\n')
    new_lines.append('.inputs ' + ' '.join(inputs) + '\n')  
    new_lines.append('.outputs ' + ' '.join(outputs) + '\n') 
    new_lines.extend(lines[start_index:])


    with open(filepath, 'w') as file:
        file.writelines(new_lines)

#新建文件夹，按照circuit描述部分的大小来决定是否更新，按照merged_1.blif,merged_2.blif这样来决定
def Incr_merge_files(folder_path):
    
    for root, dirs, files in os.walk(folder_path):
        blif_files = [file for file in files if file.endswith(".blif")]
        if blif_files:
            
            inputs = set()
            outputs = set()
            circuit_description = []
            inputs_info = []
            outputs_info = []
            circuit_info = []
            fileitem = 0

            for blif_file in blif_files:
                
                blif_file_path = os.path.join(root, blif_file)
                with open(blif_file_path, "r") as blif:
                    if len(circuit_description) >= 1000:
                        inputs_info.append(inputs)
                        outputs_info.append(outputs)
                        circuit_info.append(circuit_description)
                        inputs.clear()
                        outputs.clear()
                        circuit_description.clear() 
                        fileitem += 1

                    for line in blif:
                        # 忽略下面信息
                        if line.startswith("#"):
                            continue
                        if line.startswith(".model"):
                            continue
                        if line.startswith(".end"):
                            continue
                        
                        if line.startswith(".inputs"):
                            # 提取输入信息
                            line_inputs = line.strip().split()[1:]
                            inputs.update(line_inputs)
                        elif line.startswith(".outputs"):
                            # 提取输出信息
                            line_outputs = line.strip().split()[1:]
                            outputs.update(line_outputs)
                        else:
                            # 保存电路描述部分
                            circuit_description.append(line)

            inputs_info.append(inputs)
            outputs_info.append(outputs)
            circuit_info.append(circuit_description)

            for i in range(0,fileitem+1):
                merged_filename = 'merged' + str(i) +'.blif'
                merged_blif_path = os.path.join(root, merged_filename)
                with open(merged_blif_path, "w") as merged_blif:

                    merged_blif.write(".model "+ merged_filename +'\n')

                    merged_blif.write(".inputs ")
                    merged_blif.write(" ".join(inputs_info[i]))
                    merged_blif.write("\n")

                    merged_blif.write(".outputs ")
                    merged_blif.write(" ".join(outputs_info[i]))
                    merged_blif.write("\n")

                    for line in circuit_info[i]:
                        merged_blif.write(line)

                    merged_blif.write(".end")
            
                print("Merged BLIF file created: " + merged_blif_path)

            #return merged_blif_path

def find_duplicate_last_ports(blif_files):
    last_ports = {}
    duplicate_ports = {}
    duplicate_last_ports = set()

    index = 0
    # 遍历所有的blif文件
    for blif_file in blif_files:
        file_content_tmp = []
        file_content = []

        with open(blif_file, 'r+') as file:
            file.seek(0)
            #每一个文件中重复出现的端口
            repeated_port = []
            # 查找重复定义的.names语句的最后一个端口
            for line in file:
                if line.startswith(".names"):
                    match = re.search(r'(\S+)\s*$', line)
                    if match:
                        last_port = match.group(1)
                        if last_port in last_ports:
                            # 发现重复定义的最后一个端口
                            duplicate_last_ports.add(last_port)
                            repeated_port.append(last_port)
                            duplicate_ports[last_port] = last_port + str(index)
                            newline =  line.replace(last_port,last_port + str(index))
                            # pattern = fr'\b{re.escape(last_port)}\b'
                            # content = re.sub(pattern, last_port + str(index), content)
                            file_content_tmp.append(newline)
                            #将这一整行存储，还要生成新的对应端口，可以使用key，value，还要存储对应的文件名
                        else:
                            last_ports[last_port] = blif_file
                            file_content_tmp.append(line)
                else:
                    file_content_tmp.append(line)                   
            
            for i, line in enumerate(file_content_tmp):
                for last_port in repeated_port:
                    match = re.search(fr'\b{re.escape(last_port)}\b', line)
                    if match:
                        line = line.replace(last_port, last_port + str(index))
                        file_content_tmp[i] = line

            file.seek(0)
            file.writelines(file_content_tmp)
        index = index + 1
            
        # for key, value in duplicate_ports.items():
        #     print(key,value)
    
    return last_ports, duplicate_last_ports

if __name__ == "__main__":
    inputBlifName = sys.argv[1]
    
    folder_path =  sys.argv[2]

    cmd_1 = 'mkdir '+ folder_path
    os.system(cmd_1)

    cmd_2 = "/home/kxzhu/partition/abc_p/abc -c 'read_blif " + sys.argv[1] + "; pif /home/kxzhu/partition/abc_p/ymc_test/lib9.dsd /home/kxzhu/partition/abc_p/"+ folder_path +"'"
    os.system(cmd_2)

    print(folder_path)

    blif_files = [os.path.join(folder_path, file) for file in os.listdir(folder_path) if file.endswith('.blif')]

    last_ports, duplicate_last_ports = find_duplicate_last_ports(blif_files)

    process_blif_files(folder_path)

    Incr_merge_files(folder_path)
    
    #cmd = "/home/kxzhu/partition/abc_p/abc -c 'read_blif " + folder_path + "merged.blif; write_verilog " + folder_path + "merged.v '"
    #os.system(cmd)