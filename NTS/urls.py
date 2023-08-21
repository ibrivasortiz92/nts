from django.urls import path
from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path("main", views.main_view, name="main"),
    
    # API ROUTES
    path("ntc", views.ntc_api, name="ntc"),
]