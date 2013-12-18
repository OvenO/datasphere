
def main():
    f = open("mult_stabl.o3072930","r")
    l = f.readlines()
    new = open("q_file.txt","w")
    for i,j in enumerate(l):
        if "q" in j:
            new.write(j+"\n")
    f.close()
    new.close()
if __name__ == "__main__":
    main()
