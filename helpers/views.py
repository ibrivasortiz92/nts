from django.shortcuts import render

def handle_not_found(request, exception):
    return render(request, '_partials/not-found.html')

def handle_server_error(request):
    return render(request, '_partials/server-error.html')