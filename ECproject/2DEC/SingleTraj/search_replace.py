import os

# fix fuction works with first argument being the line to fix the second argument being the string
# you want changed and the third being the string to fix it to. it retuns the full line as it should be
# printed.
def fix(line,wrong,right):
    # find the location of the wrong string
    loc = line.find(wrong)
    # get the first corect part of the string
    first = line[:loc]
    # second part
    second = line[(loc+len(wrong)):]
    # put the right string in the sandwhic
    return first+right+second
    

def main():
    
    # string to replace
    wrong_str = '1DEC'
    # new string
    right_str = '2DEC'

    # files to omit in search and replace
    omit_files = ['.','search_replace.py']

    everything = os.listdir('.')
    to_view = []

    # get just the python or qsub script fiels in an array
    for i,j in enumerate(everything):
        if j in omit_files:
            continue

        to_view.append(j)
    
    for i,j in enumerate(to_view):
        # get the lines in a read only file and then close that file
        cur_file = open(j,'r')
        cur_lines = cur_file.readlines()
        cur_file.close()

        write_file = open(j,'w')
        for a,b in enumerate(cur_lines):
            cur_line =  b
            if wrong_str in b:
                # fix fuction works with first argument being the line to fix the second argument being the string
                # you want changed and the third being the string to fix it to. it retuns the full line as it should be
                # printed.
                cur_line = fix(b,wrong_str,right_str)
                

            write_file.write(cur_line)

        write_file.close()

if __name__ == '__main__':
    main()
