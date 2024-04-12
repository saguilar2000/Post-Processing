import re

def find_get_mass():
    with open("krome_all.f90","r") as krome_all:
        krome_all = krome_all.read()
    header = "function get_mass()"
    footer = "end function get_mass"
    start = krome_all.find(header)
    end = krome_all.find(footer)
    function = krome_all[start:end]

    return function

def build_info_dat():
    function = find_get_mass()
    pattern = re.compile(r'[Z0-9]+')
    with open("info.dat","w") as info_dat:
        for line in function.split("\n"):
            if "!" in line:
                specie = line.split("!")[1].strip()
                index = ''.join(pattern.findall(line.split("=")[0]))
                mass = line.split("=")[1].split("!")[0].strip().replace("d","e")
                info_dat.write(f"{specie}\t{index}\t{mass}\n")

def main():
    build_info_dat()

if __name__ == "__main__":
    main()