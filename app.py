import streamlit as stream
from Bio import SeqIO
import neatbio.sequtils as utils
from collections import Counter
from Bio.Seq import Seq
import requests
import io
import matplotlib.pyplot as plt
import matplotlib
stream.set_option('deprecation.showPyplotGlobalUse', False)
import numpy as np
import base64


## Sample API URL
url = "https://api.23andme.com/3/marker/?id=rs10195681,i4001358"

def main():

    @stream.cache(allow_output_mutation=True)
    def get_base64_of_bin_file(bin_file):
        with open(bin_file, 'rb') as f:
            data = f.read()
        return base64.b64encode(data).decode()

    def set_png_as_page_bg(png_file):
        bin_str = get_base64_of_bin_file(png_file)
        page_bg_img = '''
        <style>
        body {
        background-image: url("data:image/png;base64,%s");
        background-size: cover;
        }
        </style>
        ''' % bin_str

        stream.markdown(page_bg_img, unsafe_allow_html=True)
        return

    set_png_as_page_bg('background.png')


    page_bg_img = '''
    <style>
    body {
    background-image: url("https://images.unsplash.com/photo-1542281286-9e0a16bb7366");
    background-size: cover;
    }
    </style>
    '''

    stream.markdown(page_bg_img, unsafe_allow_html=True)

    stream.title("iGENOMICS")
    options = ["Introduction", "DNA Analysis","Sample API"]
    selection = stream.sidebar.radio("Choose Your Options",options)

    if selection == "Introduction":
        stream.header("A brief intro to DNA Sequence Analysis! ")

    elif selection == "DNA Analysis":
        stream.header("DNA sequence Analysis:")
        sequence=stream.file_uploader("Upload a .FASTA file for Genome analysis", type = ["fasta","fa"])

        if sequence is not None:
            data =sequence.read()
            decdata = data.decode('UTF-8')
            dnarecord= SeqIO.read(io.StringIO(decdata),"fasta")
            dnasequence= dnarecord.seq

            #print DNA reord
            details= stream.radio ("Details of the DNA by NCBI",("DNA Description", "Sequence"))
            if details == "DNA Description":
                stream.write(dnarecord.description)
            elif details == "Sequence":
                stream.write (dnarecord.seq)

            #Nucleotide
            stream.subheader ("Nucleotide Frequency :")
            dnafreq=Counter(dnasequence)
            stream.write(dnafreq)

            if stream.button("Plot frequency"):
                fig, ax = plt.subplots()
                ax.bar(dnafreq.keys(),dnafreq.values())
                stream.pyplot(fig)


            stream.subheader("DNA Composition details:")

            gc= utils.gc_content(str(dnasequence))
            at= utils.at_content(str(dnasequence))
            stream.json({"GC Content(for heat stability)": gc,"AT Content":at })

            stream.write(dnasequence)
            #Protein synthesis
            stream.subheader("Protein Synthesis:")
            ps = dnasequence.translate()
            aafreq= Counter(str(ps))

            stream.subheader("Transcription:")
            stream.write(dnasequence.transcribe())

            stream.subheader("Translation:")
            stream.write(dnasequence.translate())

            stream.subheader("Complement:")
            stream.write(dnasequence.complement())

            stream.subheader("Amino Acid frequency ::")
            stream.write(aafreq)

            # if stream.checkbox("Transcription:"):
            #     stream.write(dnasequence.transcribe())
            # elif stream.checkbox("Translation:"):
            #     stream.write(dnasequence.translate())
            # elif stream.checkbox("Complement:"):
            #     stream.write(dnasequence.complement())
            # elif stream.checkbox("Amino Acid frequency :"):
            #     stream.write(aafreq)
            #
            # elif stream.checkbox("Plot the Amino Acid frequency:"):
            #     aacolor=stream.color_picker("Pick the Amino acid color:")
            #     plt.bar(aafreq.keys(),aafreq.values(),color=aacolor)
            #     stream.pyplot()
            #
            # elif stream.checkbox("The complete Amino acid name is given as"):
            #     aaname= str(ps).replace("*","")
            #     aa3= utils.convert_1to3(aaname)
            #     stream.write(aaname)
            #     stream.write("========================")
            #     stream.write(aa3)
            #
            #
            #
            #     stream.write("========================")
            #     stream.write(utils.get_acid_name(aa3))

    elif selection == "Sample API":
        stream.header("Here is the sample API:")
        result = requests.get(url, data = {}, headers = {"Authorization": "Bearer 6b2c666cfc94ed707d6f4a4d27dcdeae"})
        stream.json({"GC Content(for heat stability)": result.text})

if __name__=='__main__':
	main()
