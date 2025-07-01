import argparse

parser = argparse.ArgumentParser(description='For scaling events by lumi')

parser.add_argument('--luminosity', action='store', default=59.0, help="Lumi scale\n DEFAULT: 59.0")
parser.add_argument('--s', action='store', default=1.0, help="signal rate\n DEFAULT: 1.0")
parser.add_argument('--b', action='store', default=1.0, help="background rate\n DEFAULT: 1.0")
parser.add_argument('--cuttype', action='store', type=str, default='a', help="cut type a or b\n DEFAULT: a")

args = parser.parse_args()
rate_list = [args.s,args.b]

# opening the file in read mode
file = open("card_modifiers/example_datacard.txt", "r")
replacement = ""
# using the for loop
for linenum, line in enumerate(file):
    #print(linenum)
    if(linenum==(17-1)):
        #print(line)
        line_mod=line
        x1 = line.split()
        x1 = x1[1:]
        #print(x1)
        for i_num, i in enumerate(x1):
            #imod=float(i)*(float(args.luminosity)/59.0)
            imod=float(rate_list[i_num])*(float(args.luminosity)/59.0)
            #print(imod)
            res = "{:.6f}".format(imod)
            #print(res)
            line_mod=line_mod.replace(i, res)
        #print(line_mod)
        replacement = replacement + line_mod
    else:
        replacement = replacement + line
            
file.close()

# opening the file in write mode
fout = open("modified_simplest_card.txt", "w")
fout.write(replacement)
fout.close()
