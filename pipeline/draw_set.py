from fpdf import FPDF
import pandas as pd
from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq

colors = [
	
	(0, 149, 239), #F3
	(60, 80, 177), #F2
	(106, 56, 179), #F1c
	(162, 36, 173), #B1c
	(243, 29, 100), #B2
	(254, 67, 60)  #B3

]


def createPDF():
	pdf = FPDF(orientation = 'P', unit = 'pt', format=(2000,500))
	pdf.set_font("Courier",'',10)
	pdf.add_page()
	return pdf

def writeln(content,r=0,g=0,b=0):
	pdf.set_text_color(r,g,b)
	pdf.cell(0,0,content)
	pdf.ln(10)

def write(content,r=0,g=0,b=0):
	pdf.set_text_color(r,g,b)
	pdf.cell(0,0,content)
	pdf.ln(0)

def bold():
	pdf.set_font("Courier","B",10)

def unbold():
	pdf.set_font("Courier",'',10)


def createB(pdf):
	bold()
	writeln("     "+" "*(info["4_position"]+len(info["4_sequence"])+1)+info["3_sequence"][::-1],*colors[3])
	writeln("     "+" "*(info["4_position"]+len(info["4_sequence"]))+"/")
	write("     "+" "*(info["4_position"])+str(Seq(info["4_sequence"]).reverse_complement()[::-1]),*colors[4])
	writeln("     "+" "*(info["4_position"])+" "*len(info["4_sequence"])+" "*(info["5_position"] - (info["4_position"]+len(info["4_sequence"])))+str(Seq(info["5_sequence"]).reverse_complement()[::-1]),*colors[5])
	unbold()
	writeln("     "+" "*(info["4_position"])+"|"*len(info["4_sequence"])+" "*(info["5_position"] - (info["4_position"]+len(info["4_sequence"])))+"|"*len(info["5_sequence"]))


def createMainSeq(pdf,sequence):
	top_seq = ""
	top_seq += "5' - "
	top_seq += str(sequence[:info["0_position"]])+" "*len(info["0_sequence"]) + str(sequence[len(info["0_sequence"])+info["0_position"]:info["1_position"]]) + " "*len(info["1_sequence"]) + str(sequence[len(info["1_sequence"])+info["1_position"] : info["3_position"]]) + " "*len(info["3_sequence"])+str(sequence[info["3_position"]+len(info["3_sequence"]) : ])

	top_seq += " - 3'"
	write(top_seq)

	bold()
	write("     "+" "*info["0_position"]+str(sequence[info["0_position"]:info["0_position"]+len(info["0_sequence"])]),*colors[0])
	write("     "+" "*info["1_position"]+str(sequence[info["1_position"]:info["1_position"]+len(info["1_sequence"])]),*colors[1])
	writeln("     "+" "*info["3_position"]+str(sequence[info["3_position"]:info["3_position"]+len(info["3_sequence"])]),*colors[3])
	unbold()


	middle_bars= "     "
	before = 0
	for i in range(6):
		middle_bars+="|"*(info[f"{i}_position"] - before)
		if i == 2 or i == 3:
			middle_bars+="|"*len(info[f"{i}_sequence"])
		else:
			middle_bars+=" "*len(info[f"{i}_sequence"])
		before = len(info[f"{i}_sequence"])+info[f"{i}_position"]
	middle_bars+="|"*(len(sequence) - before)
	writeln(middle_bars)

	bot_seq = "3' - "
	bot_seq += str(sequence.complement()[:info["2_position"]])+" "*len(info["2_sequence"])+str(sequence.complement()[len(info["2_sequence"])+info["2_position"]:info["4_position"]]) + " "*len(info["4_sequence"]) + str(sequence.complement()[len(info["4_sequence"])+info["4_position"]:info["5_position"]]) + " "*len(info["5_sequence"]) + str(sequence.complement()[len(info["5_sequence"])+info["5_position"]:])
	bot_seq += " - 5'"
	write(bot_seq)

	bold()
	write("     "+" "*info["2_position"]+str(sequence.complement()[info["2_position"]:info["2_position"]+len(info["2_sequence"])]),*colors[2])
	write("     "+" "*info["4_position"]+str(sequence.complement()[info["4_position"]:info["4_position"]+len(info["4_sequence"])]),*colors[4])
	writeln("     "+" "*info["5_position"]+str(sequence.complement()[info["5_position"]:info["5_position"]+len(info["5_sequence"])]),*colors[5])
	unbold()

def createF(pdf):
	writeln("     "+" "*info["0_position"]+"|"*len(info["0_sequence"])+" "*(info["1_position"] - (len(info["0_sequence"])+info["0_position"]))+"|"*len(info["1_sequence"]))
	
	bold()
	write("     "+" "*info["0_position"]+info["0_sequence"],*colors[0])
	writeln("     "+" "*info["1_position"]+info["1_sequence"],*colors[1])
	writeln("     "+" "*(info["1_position"]-1)+"/")
	writeln("     "+" "*(info["1_position"]-1-len(info["2_sequence"]))+str(Seq(info["2_sequence"]).reverse_complement()),*colors[2])
	unbold()



info = pd.read_csv(snakemake.input[0]).iloc[0]
reference = next(SeqIO.parse(open(snakemake.input[1],"r"),"fasta"))


pdf = createPDF()

createB(pdf)
createMainSeq(pdf,reference.seq)
createF(pdf)
#witeln(pdf,b_part+"\n"+main_part)
#pdf.multi_cell(0,10,txt="Hello")


pdf.output(snakemake.output[0])