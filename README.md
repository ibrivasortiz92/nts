## Neutron Transport Simulator (NTS)

### Video Demo: https://youtu.be/AtY52V7GDaI

### Description:

Neutron Transport Simulator (NTS) is a web application designed to provide an easy and simples way to perform neutral particle transport calculations.


### Configuration

Create virtual environment and activate it:

```bash
python3 -m venv venv
source venv/bin/activate
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Set proper environment variables (`core/.env`):

```bash
SECRET_KEY=django-secret

ALLOWED_HOST=localhost
```

### Usage

Run migrations:

```bash
python manage.py makemigrations
python manage.py migrate
```

Create admin user

```bash
python manage.py createsuperuser
```

Run Server:

```bash
python manage.py runserver
```