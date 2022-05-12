# coding:utf-8
username="XXX"
password="XXX"

import requests

auth_params = {"email": username,"password": password}
api_host = "https://www.disgenet.org/api"

## Get API key
user_key = None
with requests.Session() as s:
    try:
        r = s.post(api_host+'/auth/', data=auth_params)
        if (r.status_code == 200):
            json_response = r.json()
            user_key = json_response.get("token")
        else:
            print(r.status_code)
            print(r.text)
    except requests.exceptions.RequestException as req_ex:
        print(req_ex)
        print("Something went wrong with the request.")

