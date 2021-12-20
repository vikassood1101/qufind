from flask import Flask, render_template , redirect, url_for, request, send_from_directory
from werkzeug.utils import secure_filename
import os, json, shutil, requests
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
import qufindu as qu
import qufindv as qv
app = Flask(__name__)

import logging

logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s %(levelname)s %(message)s',
    filename= os.path.join(app.static_folder, 'logger', 'qufind.log'),
    filemode='a')

@app.route('/')
def index():
    auto_delete_file_after(getOutputFolderPath())
    return render_template("index.html")
     
@app.route('/home/')
def home():
    return redirect(url_for('index'))

@app.route('/help/')
def help():
    return render_template('help.html')

@app.route('/qufindU')
def qufindU():
    return render_template('qufindu.html')

@app.route('/submitu', methods=['GET', 'POST'])
def submitu():
    if request.method == 'POST':
        database = request.form.get('database', 'ensembl').strip()
        cpg_island = request.form.get('cpg_island', 'no').strip()
        model_type = request.form.get('model_type', 'nonoverlapping').strip()
        nucleobase = request.form.get('nucleobase', 'G').strip()
        min_stem_size = request.form.get('min_stem_size', '2').strip()
        max_stem_size = request.form.get('max_stem_size', '4').strip()
        min_loop_length = request.form.get('min_loop_length', '1').strip()
        max_loop_length = request.form.get('max_loop_length', '7').strip()
        strand = request.form.get('strand', 'positive').strip()
        bulges = request.form.get('bulges', '0').strip()
        mismatches = request.form.get('mismatches', '0').strip()

        if database in getDatabasesList()  and cpg_island in getAnswerOptions() and model_type in getOverlappingMethods() \
            and nucleobase in getNucleobase() and min_stem_size in getMinStemSize() and max_stem_size in getMaxStemSize() \
            and min_loop_length in getMinLoopLength() and max_loop_length in getMaxLoopLength() and strand in getStrands() \
            and bulges in getBulges() and mismatches in getMismatches() and not(bulges != '0' and mismatches != '0'):
            try:
                min_loop_length = int(min_loop_length)
                max_loop_length = int(max_loop_length)
                min_stem_size =  int(min_stem_size)
                max_stem_size = int(max_stem_size)
                bulges = int(bulges)
                mismatches = int(mismatches)
                random_folder_name = generateRandomString()
                valid_data_dict = {}
                invalid_data_dict = {}
                positive_dict = {}
                negative_dict = {}
                if database == 'ncbi':
                    ncbi_accession_id = request.form.get('ncbi_accession_id', '').strip()
                    if ncbi_accession_id:
                        seq_from_api = fetchDataFromAPI_id('ncbi', ncbi_accession_id)
                        if seq_from_api:
                            valid_data_dict[ ncbi_accession_id ] = seq_from_api
                else:
                    textarea_data = request.form.get('ensembl_id_textarea', '').strip()
                    if textarea_data:
                        textareDataSet = textareaDataParse(textarea_data)
                        for ids in textareDataSet:
                            seq_from_api = fetchDataFromAPI_id('ensembl', ids)
                            if seq_from_api:
                                valid_data_dict[ids] = seq_from_api
                            else:
                                invalid_data_dict[ids] = ids
                if valid_data_dict:
                    cpg_island_flag = False
                    cpg_model_value = 'nonoverlapping'
                    positive_strand_flag = True
                    negative_strand_flag =  True
                    if strand == 'positive':
                        negative_strand_flag = False
                        for key, value in valid_data_dict.items():
                            positive_dict[key] = value.upper()
                    elif strand == 'negative':
                        positive_strand_flag = False
                        for key, value in valid_data_dict.items():
                            negative_dict[key] = str(Seq(value).complement()).upper()
                    else:
                        for key, value in valid_data_dict.items():
                            positive_dict[key] = value.upper()
                            negative_dict[key] = str(Seq(value).complement()).upper()
                        
                    if cpg_island == 'yes':
                        cpg_model = request.form.get('cpg_model', 'nonoverlapping').strip()
                        cpg_island_flag = True
                        cpg_model_value = cpg_model
                    
                    if positive_dict or negative_dict:
                        if not os.path.exists(os.path.join(getOutputFolderPath(random_folder_name))):
                            os.mkdir(os.path.join(getOutputFolderPath(random_folder_name)))
                        qu.qufindu(getOutputFolderPath(random_folder_name), positive_dict, negative_dict, model_type, nucleobase, min_stem_size, max_stem_size, min_loop_length, max_loop_length, positive_strand_flag, negative_strand_flag , bulges, mismatches, cpg_island_flag, cpg_model_value)
                if invalid_data_dict:
                    if not os.path.exists(os.path.join(getOutputFolderPath(random_folder_name))):
                        os.mkdir(os.path.join(getOutputFolderPath(random_folder_name)))
                    file_name = os.path.join(getOutputFolderPath(random_folder_name), 'invalid_data.json')
                    with open(file_name, 'w') as fp_write:
                        json.dump(invalid_data_dict, fp_write)
                return redirect(url_for('viewResult', foldername = random_folder_name))
            except Exception:
                logging.exception("message")
        return render_template('error.html', msg="Something is wrong with your data!!")
    return redirect(url_for('qufindu'))

@app.route('/submitv', methods=['GET', 'POST'])
def submitv():
    if request.method == 'POST':
        model_type = request.form.get('model_type', 'nonoverlapping').strip()
        nucleobase = request.form.get('nucleobase', 'G').strip()
        min_stem_size = request.form.get('min_stem_size', '2').strip()
        max_stem_size = request.form.get('max_stem_size', '4').strip()
        min_loop_length = request.form.get('min_loop_length', '1').strip()
        max_loop_length = request.form.get('max_loop_length', '7').strip()
        bulges = request.form.get('bulges', '0').strip()
        mismatches = request.form.get('mismatches', '0').strip()

        folder_name = generateRandomString()
        valid_data_dict_1 = {}
        textarea_data = request.form.get('fasta_seq_textarea_1', '').strip()
        uploaded_data_valid_flag = False
        if textarea_data:
            response_dict = processTextarea_to_file(textarea_data)
            if 'error' in response_dict:
                return render_template('error.html', msg = response_dict['error'])
            else:
                if len(response_dict['success']) > 2:
                    return render_template('error.html', msg = 'Single sequence fasta is allowed!!')
                else:
                    valid_data_dict_1 = response_dict['success']
                    uploaded_data_valid_flag = True
        else:
            #to change
            response_dict = uploadfileFastaFormat(request, 'fasta_seq_file_1', folder_name, 'input_seq_1.txt')
            if 'error' in response_dict:
                return render_template('error.html', msg = response_dict['error'])
            else:
                valid_data_dict_1 = response_dict['success']
                uploaded_data_valid_flag = True
        
        valid_data_dict_2 = {}
        if uploaded_data_valid_flag:
            uploaded_data_valid_flag = False
            textarea_data = request.form.get('fasta_seq_textarea_2', '').strip()
            if textarea_data:
                response_dict = processTextarea_to_file(textarea_data)
                if 'error' in response_dict:
                    return render_template('error.html', msg = response_dict['error'])
                else:
                    if len(response_dict['success']) > 1:
                        return render_template('error.html', msg = 'Single sequence fasta is allowed!!')
                    else:
                        valid_data_dict_2 = response_dict['success']
                        uploaded_data_valid_flag = True
            else:
                response_dict = uploadfileFastaFormat(request, 'fasta_seq_file_2', folder_name, 'input_seq_2.txt')
                if 'error' in response_dict:
                    return render_template('error.html', msg = response_dict['error'])
                else:
                    valid_data_dict_2 = response_dict['success']
                    uploaded_data_valid_flag = True

        if valid_data_dict_1 and valid_data_dict_2 and uploaded_data_valid_flag and model_type in getOverlappingMethods() \
            and nucleobase in getNucleobase() and min_stem_size in getMinStemSize() and max_stem_size in getMaxStemSize() \
            and min_loop_length in getMinLoopLength() and max_loop_length in getMaxLoopLength() \
            and bulges in getBulges() and mismatches in getMismatches() and not(bulges != '0' and mismatches != '0'):
            try:
                min_loop_length = int(min_loop_length)
                max_loop_length = int(max_loop_length)
                min_stem_size =  int(min_stem_size)
                max_stem_size = int(max_stem_size)
                bulges = int(bulges)
                mismatches = int(mismatches)
                random_folder_name = generateRandomString()
                seq_1_dict = {}
                seq_2_dict = {}
                for key, value in valid_data_dict_1.items():
                    seq_1_dict[key] = value.upper()
                for key, value in valid_data_dict_2.items():
                    seq_2_dict[key] = value.upper()
                if seq_1_dict and seq_2_dict:
                    if not os.path.exists(os.path.join(getOutputFolderPath(random_folder_name))):
                        os.mkdir(os.path.join(getOutputFolderPath(random_folder_name)))
                    qv.qufindv(getOutputFolderPath(random_folder_name), seq_1_dict, seq_2_dict, model_type, nucleobase, min_stem_size, max_stem_size, min_loop_length, max_loop_length, bulges, mismatches)
                    return redirect(url_for('viewResultV', foldername = random_folder_name))
            except Exception:
                logging.exception("message")
        return render_template('error.html', msg="Something is wrong with your data!!")
    return redirect(url_for('qufindv'))

@app.route('/qufindV')
def qufindV(): 
    return render_template('qufindv.html')

@app.route('/qufind')
def qufind(): 
    return render_template('qufind.html')

@app.route('/submit', methods=['GET', 'POST'])
def submit():
    if request.method == 'POST':
        cpg_island = request.form.get('cpg_island', 'no').strip()
        model_type = request.form.get('model_type', 'nonoverlapping').strip()
        nucleobase = request.form.get('nucleobase', 'G').strip()
        min_stem_size = request.form.get('min_stem_size', '2').strip()
        max_stem_size = request.form.get('max_stem_size', '4').strip()
        min_loop_length = request.form.get('min_loop_length', '1').strip()
        max_loop_length = request.form.get('max_loop_length', '7').strip()
        strand = request.form.get('strand', 'positive').strip()
        bulges = request.form.get('bulges', '0').strip()
        mismatches = request.form.get('mismatches', '0').strip()

        valid_data_dict = {}
        textarea_data = request.form.get('fasta_seq_textarea', '').strip()
        uploaded_data_valid_flag = False
        if textarea_data:
            response_dict = processTextarea_to_file(textarea_data)
            if 'error' in response_dict:
                return render_template('error.html', msg = response_dict['error'])
            else:
                valid_data_dict = response_dict['success']
                uploaded_data_valid_flag = True
        else:
            response_dict = uploadfileQufind(request, generateRandomString())
            if 'error' in response_dict:
                return render_template('error.html', msg = response_dict['error'])
            else:
                valid_data_dict = response_dict['success']
                uploaded_data_valid_flag = True
            
        if valid_data_dict and uploaded_data_valid_flag and cpg_island in getAnswerOptions() and model_type in getOverlappingMethods() \
            and nucleobase in getNucleobase() and min_stem_size in getMinStemSize() and max_stem_size in getMaxStemSize() \
            and min_loop_length in getMinLoopLength() and max_loop_length in getMaxLoopLength() and strand in getStrands() \
            and bulges in getBulges() and mismatches in getMismatches() and not(bulges != '0' and mismatches != '0'):
            try:
                min_loop_length = int(min_loop_length)
                max_loop_length = int(max_loop_length)
                min_stem_size =  int(min_stem_size)
                max_stem_size = int(max_stem_size)
                bulges = int(bulges)
                mismatches = int(mismatches)
                random_folder_name = generateRandomString()
                positive_dict = {}
                negative_dict = {}
                cpg_island_flag = False
                cpg_model_value = 'nonoverlapping'
                positive_strand_flag = True
                negative_strand_flag =  True
                if strand == 'positive':
                    negative_strand_flag = False
                    for key, value in valid_data_dict.items():
                        positive_dict[key] = value.upper()
                elif strand == 'negative':
                    positive_strand_flag = False
                    for key, value in valid_data_dict.items():
                        negative_dict[key] = str(Seq(value).complement()).upper()
                else:
                    for key, value in valid_data_dict.items():
                        positive_dict[key] = value.upper()
                        negative_dict[key] = str(Seq(value).complement()).upper()
                    
                if cpg_island == 'yes':
                    cpg_model = request.form.get('cpg_model', 'nonoverlapping').strip()
                    cpg_island_flag = True
                    cpg_model_value = cpg_model
                
                if positive_dict or negative_dict:
                    if not os.path.exists(os.path.join(getOutputFolderPath(random_folder_name))):
                        os.mkdir(os.path.join(getOutputFolderPath(random_folder_name)))
                    qu.qufindu(getOutputFolderPath(random_folder_name), positive_dict, negative_dict, model_type, nucleobase, min_stem_size, max_stem_size, min_loop_length, max_loop_length, positive_strand_flag, negative_strand_flag , bulges, mismatches, cpg_island_flag, cpg_model_value)
            
                return redirect(url_for('viewResult', foldername = random_folder_name))
            except Exception:
                logging.exception("message")
        return render_template('error.html', msg="Something is wrong with your data!!")
    return redirect(url_for('qufind'))


def uploadfileFastaFormat(request, uploaded_file_name, folder_name, input_file_name):
    Upload_folder = getOutputFolderPath(folder_name)
    try:
        data_dict = {}
        errorMsgList = []
        fileUploadResponseDict = {}
        file = request.files[uploaded_file_name]
        file_upload_error = validateUploadedFile(file)
        if file_upload_error:
            errorMsgList.append('Sequence File: ' + file_upload_error)
        else:
            if not os.path.exists(Upload_folder):
                os.mkdir(Upload_folder)
            filename = secure_filename(input_file_name)
            filePath = os.path.join(Upload_folder, filename)
            file.save(filePath)
            if os.stat(filePath).st_size > 0:
                try:
                    #Single Sequence Fasta
                    fasta_file_data_obj = SeqIO.read(open(filePath), "fasta")
                    data_dict[fasta_file_data_obj.id] = str(fasta_file_data_obj.seq)
                except:
                    errorMsgList.append('Sequence File is not in proper fasta format.')
            else:
                errorMsgList.append('Sequence File: File is Empty!!')
        
        if errorMsgList and not data_dict:
            fileUploadResponseDict['error'] = '</br>'.join(errorMsgList)
        else:
            fileUploadResponseDict['success'] = data_dict
    except Exception:
        logging.exception('message')
        fileUploadResponseDict['error'] = 'Sequence File can not be processed.'
    if os.path.exists(Upload_folder):
        shutil.rmtree(Upload_folder)
    return fileUploadResponseDict

def uploadfileQufind(request, folder_name):
    Upload_folder = getOutputFolderPath(folder_name)
    try:
        data_dict = {}
        errorMsgList = []
        fileUploadResponseDict = {}
        file = request.files['qufind_fasta_seq_file']
        file_upload_error = validateUploadedFile(file)
        if file_upload_error:
            errorMsgList.append('Sequence File: ' + file_upload_error)
        else:
            if not os.path.exists(Upload_folder):
                os.mkdir(Upload_folder)
            filename = 'input_fasta.txt'
            filePath = os.path.join(Upload_folder, filename)
            file.save(filePath)
            if os.stat(filePath).st_size > 0:
                try:
                    fasta_file_data_obj = SeqIO.parse(filePath, "fasta")
                    for seq_record in fasta_file_data_obj:
                        data_dict[seq_record.id] = str(seq_record.seq)
                except:
                    errorMsgList.append('Sequence File is not in proper Fasta format.')
            else:
                errorMsgList.append('Sequence File: File is Empty!!')
        
        if errorMsgList and not data_dict:
            fileUploadResponseDict['error'] = '</br>'.join(errorMsgList)
        else:
            fileUploadResponseDict['success'] = data_dict
    except Exception:
        logging.exception('message')
        fileUploadResponseDict['error'] = 'Sequence File can not be processed.'
    if os.path.exists(Upload_folder):
        shutil.rmtree(Upload_folder)
    return fileUploadResponseDict

def validateUploadedFile(file_obj):
    error_msg = ''
    if file_obj.filename == '':
        error_msg = 'No File Selected!!'
    elif file_obj and file_obj.filename:
        if not allowed_file(file_obj.filename):
            error_msg = 'File format not allowed!!'
    return error_msg

def allowed_file(filename):
    ALLOWED_EXTENSIONS = { 'txt', 'fa', 'fasta', 'fsa', 'fna', 'faa', 'ffn', 'frn'}
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/contactus/')
def contactus():
    return render_template('contactus.html')

def getDatabasesList():
    return ['ensembl', 'ncbi']

def getAnswerOptions():
    return ['yes', 'no']

def getOverlappingMethods():
    return ['overlapping', 'nonoverlapping']

def getNucleobase():
    return ['G', 'C', 'A', 'T']

def getMinStemSize():
    return ['2' , '3']

def getMaxStemSize():
    return ['4' , '5']

def getMinLoopLength():
    return [ str(i) for i in range(1, 11)]

def getMaxLoopLength():
    return [ str(i) for i in range(5, 31)]

def getStrands():
    return ['positive', 'negative', 'both']

def getBulges():
    return [ str(i) for i in range(0, 6)]

def getMismatches():
    return ['0', '1', '2']

def read_fasta_file(fastaFilePath):
    positive_dict = {}
    negative_dict = {}
    for seq_record in SeqIO.parse(fastaFilePath, "fasta"):
        positive_dict[seq_record.id] = str(seq_record.seq).upper()
        negative_dict[seq_record.id] = str(seq_record.seq.complement()).upper()
    return positive_dict, negative_dict

@app.route('/viewResultV/<foldername>')
def viewResultV(foldername = None):
    #auto_delete_file_after(getOutputFolderPath(''))
    foldername = secure_filename(foldername)
    result_folder = getOutputFolderPath(foldername)
    json_file_dict = {}
    plot_image_encoded = ''
    if os.path.exists(result_folder):
        for file_ in os.listdir(result_folder):
            current_file_path = os.path.join(result_folder, file_)
            if os.path.isfile(current_file_path) and (file_.endswith('unique_data.json') or file_.endswith('common_data.json')):
                json_data_dict = processJsonFile(current_file_path)
                if json_data_dict:
                    del json_data_dict['header']
                    inner_data_dict = {}
                    for seq_id, table_data in json_data_dict.items():
                        if table_data:
                            inner_data_list = []
                            inner_data_list.append(
                                [
                                    "G4 sequence",
                                    "Length",
                                    "Start Position",
                                    "End Position"
                                ]
                            )
                            for key, value in table_data.items():
                                for ele in value:
                                    start_position = int(key)
                                    end_position = start_position + len(ele)
                                    inner_data_list.append(
                                        [ ele, len(ele), start_position, end_position ]
                                    )
                            inner_data_dict[seq_id] = inner_data_list
                    if inner_data_dict:
                        if 'unique' in file_:
                            data_dict = json_file_dict.get('unique_quadruplexes', {})
                            data_dict.update(inner_data_dict)
                            json_file_dict[ 'unique_quadruplexes' ] = data_dict
                        else:
                            data_dict = json_file_dict.get('common_quadruplexes', {})
                            data_dict.update(inner_data_dict)
                            json_file_dict[ 'common_quadruplexes' ] = data_dict

        plot_image_path = os.path.join(result_folder, 'plot_over.png')
        if os.path.exists(plot_image_path):
            plot_image_encoded = encodeImageToBase64(plot_image_path)
            
        if json_file_dict:
            return render_template('outputv.html', plot_image = plot_image_encoded, foldername = foldername, folder_data_dict = json_file_dict)
        return render_template('error.html', msg ='Sorry!! No Quadruplexes were found in your sequence.')        
    return render_template('error.html', msg ='Sorry-- Data not found/has been erased!!!')

@app.route('/viewResult/<foldername>')
def viewResult(foldername = None):
    #auto_delete_file_after(getOutputFolderPath(''))
    foldername = secure_filename(foldername)
    result_folder = getOutputFolderPath(foldername)
    json_file_dict = {}
    if os.path.exists(result_folder):
        resultGenerated = False
        for file_ in os.listdir(result_folder):
            current_file_path = os.path.join(result_folder, file_)
            if os.path.isfile(current_file_path) and file_.endswith('json') and ('cpg' in file_ or 'motif' in file_):
                json_data_dict = processJsonFile(current_file_path)
                if json_data_dict:
                    resultGenerated = True
                    del json_data_dict['header']
                    inner_data_dict = {}
                    if 'cpg' in file_:
                        for seq_id, table_data in json_data_dict.items():
                            inner_data_list = []
                            inner_data_list.append(
                                [
                                    "Region Detected(Positions)",
                                    "GC Content",
                                    "CpG Ratio"
                                ]
                            )
                            for key, value in table_data.items():
                                inner_data_list.append(
                                    [ key + '-' + str(value[0]), value[1], value[2]]
                                )
                            if 'positive' in file_:
                                inner_data_dict['positive_cpg'] = inner_data_list
                            else:
                                inner_data_dict['negative_cpg'] = inner_data_list
                            data_dict = json_file_dict.get(seq_id, {})
                            data_dict.update(inner_data_dict)
                            json_file_dict[ seq_id ] = data_dict
                            current_png_path = os.path.join(result_folder, seq_id + '.png')
                            if os.path.exists(current_png_path) and 'png_data' not in json_file_dict[ seq_id ]:
                                encoded_png = encodeImageToBase64(current_png_path)
                                if encoded_png:
                                    json_file_dict[ seq_id ].update({'png_data' : encoded_png})
                    else:
                        for seq_id, table_data in json_data_dict.items():
                            inner_data_list = []
                            inner_data_list.append(
                                [
                                    "G4 sequence",
                                    "Length",
                                    "Start Position",
                                    "End Position"
                                ]
                            )
                            for key, value in table_data.items():
                                for ele in value:
                                    start_position = int(key)
                                    end_position = start_position + len(ele)
                                    inner_data_list.append(
                                        [ ele, len(ele), start_position, end_position ]
                                    )
                            if 'positive' in file_:
                                inner_data_dict['positive_motif'] = inner_data_list
                            else:
                                inner_data_dict['negative_motif'] = inner_data_list
                            data_dict = json_file_dict.get(seq_id, {})
                            data_dict.update(inner_data_dict)
                            json_file_dict[ seq_id ] = data_dict
                            current_png_path = os.path.join(result_folder, seq_id + '.png')
                            if os.path.exists(current_png_path) and 'png_data' not in json_file_dict[ seq_id ]:
                                encoded_png = encodeImageToBase64(current_png_path)
                                if encoded_png:
                                    json_file_dict[ seq_id ].update({'png_data' : encoded_png})
        if resultGenerated:
            if json_file_dict:
                return render_template('output.html', foldername = foldername, folder_data_dict = json_file_dict)
        return render_template('error.html', msg ='Sorry!! No Quadruplexes were found in your sequence.')        
    return render_template('error.html', msg ='Sorry-- Data not found/has been erased!!!')

def processJsonFile(filePath):
    data_dict = {}
    with open(filePath) as file_read:
        data_dict = json.load(file_read)
    return data_dict

@app.template_filter('tableHeader')
def tableOutputHeaderNames(input_str):
    header_dict = {
        'negative_cpg' : 'CPG Negative',
        'negative_motif' : 'Motif Negative',
        'positive_cpg' : 'CPG Positive',
        'positive_motif' : 'Motif Positive',
        'seq_1_motif': '1st Sequence',
        'seq_2_motif': '2nd Sequence',
        'unique_quadruplexes': 'Unique Quadruplexes',
        'common_quadruplexes': 'Common Quadruplexes'
    }
    if input_str not in header_dict:
        return input_str
    return header_dict[input_str]


@app.route('/downloadFile', methods=['GET', 'POST'])
def downloadFile():
    if request.method == 'POST':
        try:
            foldername = request.form.get('foldername', '').strip()
            if foldername:
                text_file_list = []
                if '#' in foldername and foldername.count('#') == 2:
                    foldername_list = foldername.split('#')
                    output_folder_path = foldername_list[0]
                    output_id = foldername_list[1]
                    output_file_type = foldername_list[2]
                    foldername = 'Output_' + tableOutputHeaderNames(output_file_type).replace(' ','_') + '.csv'
                    id_folder_file = os.path.join(getOutputFolderPath(output_folder_path), f'output_{output_file_type}.json')
                    if os.path.exists(id_folder_file):
                        json_data_dict = processJsonFile(id_folder_file)
                        if json_data_dict:
                            text_file_list.append(','.join(json_data_dict['header']) + '\n')
                            del json_data_dict['header']
                            for key, value in json_data_dict.items():
                                if key == output_id:
                                    foldername = key + '_' + foldername
                                    if 'cpg' in output_file_type:
                                        for i_key, i_value in value.items():
                                            text_file_list.append(f'{key},{i_key}-{i_value[0]},{i_value[1]},{i_value[2]}\n')
                                    else:
                                        for i_key, i_value in value.items():
                                            for ele in i_value:
                                                text_file_list.append(f'{key},{ele},{len(ele)},{i_key},{int(i_key) + len(ele)}\n')
                else:
                    program_result_folder = os.path.join(getOutputFolderPath(foldername))
                    foldername = 'Combined_Output.txt'
                    if os.path.exists(program_result_folder):
                        for output_files in os.listdir(program_result_folder):
                            if output_files.endswith('.json'):
                                output_file_name = output_files.replace('output_','').replace('.json', '')
                                text_file_list.append(f'\n\t\t\tOutput {tableOutputHeaderNames(output_file_name)}\n\n')
                                json_file = os.path.join(program_result_folder, output_files)
                                json_data_dict = processJsonFile(json_file)
                                if json_data_dict:
                                    text_file_list.append('\t'.join(json_data_dict['header']) + '\n')
                                    del json_data_dict['header']
                                    for key, value in json_data_dict.items():
                                        if 'cpg' in output_files:
                                            for i_key, i_value in value.items():
                                                text_file_list.append(f'{key}\t{i_key}-{i_value[0]}\t{i_value[1]}\t{i_value[2]}\n')
                                        else:
                                            for i_key, i_value in value.items():
                                                for ele in i_value:
                                                    text_file_list.append(f'{key}\t{ele}\t{len(ele)}\t{i_key}\t{int(i_key) + len(ele)}\n')
                        
                if text_file_list:        
                    r = app.response_class(text_file_list, mimetype='text/plain')
                    r.headers.set('Content-Disposition', 'attachment', filename = foldername)
                    return r
        except Exception:
            logging.exception("message")
        return render_template('error.html', msg = 'Result Download Error!!!!')
    return redirect(url_for('index'))

@app.route('/error', methods = ['GET', 'POST'])
def errorPage():
    errmsg = ''
    if request.method == 'POST':
        errmsg = request.form['errorInputJS']
    else:
        errmsg = request.args.get('msg', '')
    if errmsg:
        return render_template('error.html', msg = errmsg)
    return render_template('error.html')

@app.route('/loadDownloadFile/<file_name>')
def loadDownloadFile(file_name = None):
    if file_name:
        folder_location = os.path.join(app.static_folder, 'examples')
        if os.path.exists(os.path.join(folder_location, file_name + '.fasta')):
            return send_from_directory(folder_location, file_name + '.fasta', as_attachment= True)
    return redirect(url_for('index'))

##### Ajax call ######
@app.route('/getExample/<module_type>')
def getExample(module_type):
    module_type = secure_filename(module_type)
    sequenceExampleDict = {}
    errorMsg = 'Error while fetching example'
    successMsg = ''
    try:
        if module_type in ['qufindu', 'qufind', 'qufindv_1', 'qufindv_2']:
            data_list = []
            example_folder_path = os.path.join(app.static_folder, 'examples')
            example_path = os.path.join(example_folder_path, module_type + '.txt')
            if module_type == 'qufindu':
                with open(example_path) as f_read:
                    for line in f_read:
                        line = line.strip()
                        data_list.append(line)
            elif 'qufindv' in module_type:
                fasta_file_data_obj = SeqIO.read(open(example_path), "fasta")
                data_list.append('>' + fasta_file_data_obj.id)
                data_list.append(str(fasta_file_data_obj.seq))
            else:
                fasta_file_data_obj = SeqIO.parse(example_path, "fasta")
                for seq_record in fasta_file_data_obj:
                    data_list.append('>' + seq_record.id)
                    data_list.append(str(seq_record.seq))
            if data_list:
                successMsg = '\r\n'.join(data_list)
            
    except Exception:
        logging.exception('message')
    if errorMsg and not successMsg:
        sequenceExampleDict['error'] = errorMsg
    else:
        sequenceExampleDict['success'] = successMsg
    return json.dumps(sequenceExampleDict)

def getOutputFolderPath(folder_name = ''):
    return os.path.join(app.static_folder, 'results', folder_name)

def fetchDataFromAPI_id(database, id):
    respons = ''
    try:
        if database == 'ncbi':
            # server = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={id}&rettype=fasta&retmode=text'
            # resp = requests.get(server)
            # if resp.ok:
            #     data_dict = resp.text
            #     print(data_dict)
            Entrez.email = 'vino@anuefa.com'
            text_file_fp = Entrez.efetch(db="nuccore", id=id, rettype="fasta", retmode="text")
            #To discard Accessnion id Header returned in file
            text_file_fp.readline()
            for line in text_file_fp:
                respons += line.strip()
            text_file_fp.close()
        else:
            server = "https://rest.ensembl.org"
            ext = "/sequence/id/" + id
            if 'ens' in id.lower(): 
                ext += '?type=genomic'
            else:
                ext += "?type=cds"
            headers = { "Content-Type" : "application/json"}
            resp = requests.get(server+ext, headers=headers)
            if resp.ok:
                data_dict = resp.json()
                respons = data_dict['seq']
            resp.close()
    except Exception:
        logging.exception('message')
    return respons

def processTextarea_to_file(textarea_data_string):
    response_dict = {}
    textarea_data_string = textarea_data_string.replace('\r\n', ' ')
    textarea_data_string = textarea_data_string.replace('\r', ' ')
    textarea_data_string = textarea_data_string.replace('\n', ' ')
    textarea_data_list = textarea_data_string.split()
    if '>' in textarea_data_list[0]:
        seq_data_dict = {}
        header = ''
        for ele in textarea_data_list:
            if '>' == ele:
                seq_data_dict = {}
                break
            elif '>' in ele:
                header = ele[1:]
                seq_data_dict[ header ] = ''
            else:
                seq_data_dict[header] += ele
        
        for key, value in seq_data_dict.items():
            if not value:
                seq_data_dict = {}
                break

        if seq_data_dict:
            response_dict['success'] = seq_data_dict
        else:
            response_dict['error'] = 'Something is wrong with your data.'
    else:
        response_dict['error'] = 'Data is not in Fasta format.'
    return response_dict

def textareaDataParse(textarea_data_string):
    textareaDataSet = set()
    textarea_data_string = textarea_data_string.replace('\r\n', ', ')
    textarea_data_string = textarea_data_string.replace('\r', ', ')
    textarea_data_string = textarea_data_string.replace('\n', ', ')
    textarea_data_string = adjust_Tabs_Spaces(textarea_data_string)
    if ',' in textarea_data_string:
        textareaDataSet = { ele.strip() for ele in textarea_data_string.split(',') if ele.strip() }
    else:                   
        textareaDataSet.add(textarea_data_string.strip())
    return textareaDataSet

def encodeImageToBase64(imageFilePath):
    encoded_image = ''
    if os.path.exists(imageFilePath):
        import base64
        fileimage_read = open(imageFilePath, 'rb')
        binary_fileimage_read = fileimage_read.read()
        encoded_image = base64.b64encode(binary_fileimage_read).decode('utf-8')
        if encoded_image:
            encoded_image = 'data:image/png;base64,'+ encoded_image
        fileimage_read.close()
    return encoded_image

def adjust_Tabs_Spaces(valueString):
    valueString = [i for j in valueString.split() for i in (j, ' ')][:-1] #To adjust spaces to get same space length
    return ''.join(valueString)

def not_allowed_file(filename):
    NOT_ALLOWED_EXTENSIONS = { 'pdf', 'jpg', 'jpeg', 'png', 'svg', 'doc', 'docx', 'xls', 'xlsx', 'exe', 'zip','rar','dmg'}
    return not( '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in NOT_ALLOWED_EXTENSIONS )

def generateRandomString(len_ = 40):
    import string, random
    lettersAndDigits = string.ascii_letters + string.digits
    return ''.join((random.choice(lettersAndDigits) for i in range(len_)))

def auto_delete_file_after(folderpath):
    import time
    time_now = time.time()
    for filename in os.listdir(folderpath):
        current_file_folder = os.path.join(folderpath, filename)
        #if os.path.getmtime(current_file_folder) < time_now - 7 * 86400:
        if os.stat(current_file_folder).st_mtime < time_now - 1 * 86400:
            shutil.rmtree(current_file_folder)

#Ajax Call
@app.route('/drawCPGPlot')
def drawCPG_Quad_Plot():
    result_json_dict = { 'error' : 'Some Error Occurred!!'}
    query_param = request.args.get('q', '').strip()
    if query_param:
        try:
            query_param_list = query_param.split('#',4)
            folder_name = secure_filename(query_param_list[0])
            strand_type = secure_filename(query_param_list[1])
            query_id = query_param_list[2]
            start_index = secure_filename(query_param_list[3])
            end_index = secure_filename(query_param_list[4])
            if all([folder_name, strand_type, query_id, start_index, end_index]):
                outputfolderPath = getOutputFolderPath(folder_name)
                motif_main_list = []
                validParams = False
                for file_ in os.listdir(outputfolderPath):
                    if 'output_' + strand_type in file_ and file_.endswith('.json'):
                        if 'motif' in file_:
                            motif_json_dict = {}
                            with open(os.path.join(outputfolderPath, file_)) as f_read:
                                motif_json_dict = json.load(f_read)
                            del motif_json_dict['header']
                            if query_id in motif_json_dict:
                                for key, value in motif_json_dict[query_id].items():
                                    for item in value:
                                        motif_main_list.append([int(key), int(key)+len(item)])
                        else:
                            cpg_json_dict = {}
                            with open(os.path.join(outputfolderPath, file_)) as f_read:
                                cpg_json_dict = json.load(f_read)
                            del cpg_json_dict['header']
                            if query_id in cpg_json_dict:
                                if start_index in cpg_json_dict[query_id]:
                                    if end_index == str(cpg_json_dict[query_id][start_index][0]):
                                        validParams = True
        
                if motif_main_list and validParams:
                    import matplotlib.pyplot as plt
                    start_index = int(start_index)
                    end_index = int(end_index)
                    g4_start = [] 
                    g4_end = []
                    import tempfile, base64
                    tmpfile = tempfile.TemporaryFile()
                    for ele2 in motif_main_list:
                        count = 0
                        if ((ele2[0] >= start_index) and (ele2[1] <= end_index)):  
                            g4_start.append(ele2[0])
                            g4_end.append(ele2[1])
                            y = range(0, len(g4_start))
                            plt.figure(figsize=(4,3))
                            plt.hlines(y, g4_start, g4_end, color = 'red', label = 'Quadruplexes')            
                            for x1,x2 in zip(g4_start, g4_end):
                                plt.text(x1 + 0.5, y[count], s = (x1,x2), horizontalalignment = 'right', fontweight = 'bold', verticalalignment='bottom',fontsize=4, color = 'blue')
                                count = count + 1
                            plt.xlim([start_index, end_index])
                            plt.yticks([])
                            plt.xticks([start_index, end_index])
                            plt.xlabel("Positions of Quadruplexes", fontsize = 6)
                            plt.legend(prop = {'style': 'italic'}, loc='upper left', bbox_to_anchor=(1, 0.5))
                            plt.savefig(tmpfile, pad_inches = 0.1,dpi = 200, quality = 95,bbox_inches = "tight")
                            plt.close()
                            tmpfile.seek(0)
                    encoded_image = base64.b64encode(tmpfile.read()).decode('utf-8')
                    tmpfile.close()
                    if encoded_image:
                        result_json_dict = {}
                        result_json_dict['img'] = 'data:image/png;base64,'+ encoded_image
                    if not g4_start and not g4_end:
                        result_json_dict = {}
                        result_json_dict['msg'] = 'No quadruplexes were found in this cpg island'
        except Exception:
            logging.exception('message')
    return json.dumps(result_json_dict)

                
if __name__ == "__main__":    
    app.run()