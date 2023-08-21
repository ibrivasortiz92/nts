from django.shortcuts import render, redirect
from django.contrib import messages
from django.contrib.auth import authenticate, login, logout
from django.urls import reverse
from helpers.decorators import auth_user_should_not_access


@auth_user_should_not_access
def login_view(request):
    if request.method == 'POST':
        user_email = request.POST.get('email')
        user_password = request.POST.get('password')

        user = authenticate(request, email=user_email, password=user_password)

        if not user:
            messages.error(request, 'Invalid credentials')
            return render(request, 'authentication/login.html', {'user_email': user_email}, status=409)

        login(request, user)
        return redirect(reverse('main'))

    return render(request, 'authentication/login.html')



def logout_view(request):
    logout(request)
    messages.success(request, 'Log out successfully')
    return redirect(reverse('login'))
