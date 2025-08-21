from django.shortcuts import render
from django.http import HttpResponse
from django.http import JsonResponse


def home(request):
    return render(request, "home.html")

def ketcher(request):
    return render(request, "index.html")

def getshifts(request):
    mol_smiles = request.GET['molecule'].split(".")[0]
    print(mol_smiles)
    data = {'molecule': mol_smiles} #placeholder
    return JsonResponse(data)
