from django.shortcuts import render
from django.contrib.auth.decorators import login_required
from django.views.decorators.csrf import csrf_protect
from django.http import JsonResponse
from django.conf import settings
import json
import subprocess
import math

# INITIAL VIEW
def index(request):
    return render(request, "NTS/index.html")


# MAIN VIEW
@login_required
def main_view(request):
    return render(request, "NTS/main.html")


# MAIN API 
@csrf_protect
def ntc_api(request):

    # CHECK POST METHOD
    if request.method != "POST":
        return JsonResponse({
            "STATUS": 8
        }, status=400)
    
    data = json.loads(request.body)
    quad = data['QUAD']['miu']

    # CALCULATION DETAILS:
    # NAME
    detail_name = data["Name"]
    # METHOD
    detail_method = ""
    if data["METH"] == "dd_method":
        detail_method = "DD"
    elif data["METH"] == "cd_method":
        detail_method = "CD"
    elif data["METH"] == "ld_method":
        detail_method = "LD"
    elif data["METH"] == "cn_method":
        detail_method = "CN"
    elif data["METH"] == "rm_cn_method":
        detail_method = "RM-CN"
    elif data["METH"] == "rm_lln_method":
        detail_method = "RM-LLN"
    # QUADRATURE
    detail_quadrature = int(math.sqrt(1 + 2 * len(data["QUAD"]["miu"])) - 1)
    # MESH
    ncx = 0
    for i in range(data["XREG"]):
        ncx = ncx + data["XDOM"][i]["nc"];
    ncy = 0
    for j in range(data["YREG"]):
        ncy = ncy + data["YDOM"][j]["nc"];
    detail_mesh = f"{ncy} x {ncx}"
    # TOLERANCE
    detail_tolerance = data["TOL"]
    detail = {
        'name': detail_name,
        'method': detail_method,
        'quadrature': detail_quadrature,
        'mesh': detail_mesh,
        'tolerance': detail_tolerance
    }

    # CONSTRUCT INPUT
    input_file = str(int(math.sqrt(1 + 2 * len(quad)) - 1)) + '\n'
    input_file = input_file + str(data['ZN']) + '\n'
    for zone in data["ZON"]:
        input_file = input_file + str(zone['st']) + ' ' + str(zone['ss']) + '\n'
    input_file = input_file + str(data['XREG']) + '\n'
    for r in data['XDOM']:
        input_file = input_file + str(r['len']) + ' ' + str(r['nc']) + '\n'
    input_file = input_file + str(data['YREG']) + '\n'
    for r in data['YDOM']:
        input_file = input_file + str(r['len']) + ' ' + str(r['nc']) + '\n'
    for j in range(data['YREG']):
        for i in range(data['XREG']):
            input_file = input_file + str(1 + data['ZMAP'][j][i]) + ' '
        input_file = input_file + '\n'
    for j in range(data['YREG']):
        for i in range(data['XREG']):
            input_file = input_file + str(data['QMAP'][j][i]) + ' '
        input_file = input_file + '\n'
    input_file = input_file + str(data['BC']['left']) + ' ' + str(data['BC']['bottom']) + ' '
    input_file = input_file + str(data['BC']['right']) + ' ' + str(data['BC']['top']) + '\n'
    input_file = input_file + str(data['TOL']) + '\n'

    # RUN METHOD
    try:
        results = False
        n_loop = 1
        cpu_time = 0
        if data['cpu_loop']: 
            n_loop = 10
        for i in range(n_loop):
            output = {}
            output["stdout"] = False

            # RUN DD METHOD
            if data["METH"] == "dd_method":
                output = subprocess.run([f"{settings.BASE_DIR}/programs/linux/NTS_DD"],input=input_file.encode(), capture_output=True)
            
            # RUN CD METHOD
            elif data["METH"] == "cd_method":
                output = subprocess.run([f"{settings.BASE_DIR}/programs/linux/NTS_CD"],input=input_file.encode(), capture_output=True)

            # RUN LD METHOD
            elif data["METH"] == "ld_method":
                output = subprocess.run([f"{settings.BASE_DIR}/programs/linux/NTS_LD"],input=input_file.encode(), capture_output=True)

            # RUN CN METHOD
            elif data["METH"] == "cn_method":
                output = subprocess.run([f"{settings.BASE_DIR}/programs/linux/CN"],input=input_file.encode(), capture_output=True)

            # RUN RM-CN METHOD
            elif data["METH"] == "rm_cn_method":
                output = subprocess.run([f"{settings.BASE_DIR}/programs/linux/NTS_RM-CN"],input=input_file.encode(), capture_output=True)

            # RUN RM-LLN METHOD
            elif data["METH"] == "rm_lln_method":
                output = subprocess.run([f"{settings.BASE_DIR}/programs/linux/NTS_RM-LLN"],input=input_file.encode(), capture_output=True)

            if output.stdout:
                results = json.loads(output.stdout.decode())
                if results["STATUS"] == 0:
                    cpu_time = cpu_time + results['CPU']
        
        # RESULTS
        if output.stdout:

            results = json.loads(output.stdout.decode())

            # CHECK RESULT STATUS
            if results["STATUS"] != 0:
                return JsonResponse({
                    "STATUS": results['STATUS']
                })

        # SEND RESULTS
        return JsonResponse({
            "STATUS": results['STATUS'],
            "ITER": results["ITER"],
            "CPU": cpu_time / n_loop,
            "MFLUX": results["MFLUX"],
            "MFLOW": results["MFLOW"],
            "XFLOW": results["XFLOW"],
            "YFLOW": results["YFLOW"]
        })

    except:

        # SOMETHING WENT WRONG WITH THE SERVER
        return JsonResponse({
            "STATUS": 7
        })


