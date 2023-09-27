import os
from pathlib import Path

def find_command_in_file(start, end, file):
    len_start = len(start)
    len_end = len(end)
    i = 0
    start_word = ""
    end_word = ""
    command = start
    found_start_word = False
    found_end_word = False
    with open(file, 'r') as fp:
        while True:
            # Read file character by character.
            char = fp.read(1)
            if not char:
                return None
            # Add to start and end words.
            start_word += char
            end_word += char
            if found_start_word == True and found_end_word == False:
                # Start adding characters to command.
                command += char
            if i >= len_start:
                # Check for start word.
                start_word = start_word[1:len_start+1]
                if found_start_word == False and start_word == start:
                    found_start_word = True
            if i >= len_end:
                # Check for end word.
                end_word = end_word[1:len_end+1]
                if found_end_word == False and end_word == end:
                    found_end_word = True
            if found_end_word == True:
                # Return command after removing erroneously printed message.
                command = command[0:len(command)-len_end]
                if start == "ifort":
                    f = open("ifort.txt", "w")
                    f.write(command)
                    f.close()
                if start == "LINK":
                    command = command.replace("   Creating library standardU.lib and object standardU.exp\n", "")
                    f = open("link.txt", "w")
                    f.write(command)
                    f.close()
                if start == "export.sym: ":
                    command = command.replace("export.sym: ", "")
                    command = "LIBRARY standardU.dll\n" + "EXPORTS\n" + command
                    f = open("export.def", "w")
                    f.write(command)
                    f.close()
                return command
            i += 1
            
        
def find_file(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)
            
def find_word_in_file(word, file):
    with open(file, 'r') as fp:
        lines = fp.readlines()
        for line in lines:
            if line.find(word) != -1:
                return line
    return None

abaqus_v6 = find_file("abaqus_v6.env", "C:/Simulia/")
print(abaqus_v6)

win86_64 = find_file("win86_64.env", "C:/Simulia/")
print(win86_64)

print("Checking user subroutines verification...")

os.system("abaqus verify -user_std > verify.txt")
check_file = Path("verify.txt")
if check_file.is_file():
    if(find_word_in_file("Verification procedure complete", "verify.txt")):
        if(find_word_in_file("PASS", "verify.txt")):
            print("Verification of std user subroutines successful...")
            user_subroutines = True
        else:
            print("Verification of std user subroutines failed...")
            user_subroutines = False
            
if user_subroutines:
    print("Making subroutine library...")
    try:
        os.remove("standardU.dll")
        os.remove("std_user-std.obj")
    except OSError:
        pass
    os.system("abaqus make library=std_user.for > make.txt")
    print("Subroutine library made...")
    check_file = Path("make.txt")
    if check_file.is_file():
        compile_command = find_command_in_file("ifort", "End Compiling Abaqus/Standard User Subroutines", "make.txt")
        link_command = find_command_in_file("LINK", "End Linking Abaqus/Standard User Subroutines", "make.txt")  
        export_command = find_command_in_file("export.sym: ", "Linking std_user.obj into user subroutine shared library standardU.dll", "make.txt")
    try:
        os.remove("standardU.dll")
        os.remove("std_user-std.obj")
    except OSError:
        pass  
    os.system(compile_command)
    os.system(link_command)
    print("Compile and link commands all work!")
