#!/usr/bin/env python
# coding: utf-8

# # Initializing the SysML v2 API

# In[2]:


from __future__ import print_function

import time
import requests
from pprint import pprint
import pandas as pd
import json
from datetime import datetime

#host = "<specify protocol://host:port of the server that is a provider of the SysML v2 REST/HTTP API"
host = "https://flexo.systemsmodeling.com/sysmlv2"


# # Get projects

# In[4]:


projects_url = f"{host}/projects" 
projects_response = requests.get(projects_url)

if projects_response.status_code == 200:
    projects = projects_response.json()
    projects_data = list(map(lambda b: {'Project Name':b['name'], 'Project ID':b['@id']}, projects))
    df = pd.DataFrame.from_records(projects_data)
    # df = pd.DataFrame.from_records(projects_data).style.hide(axis='index')
    pprint(df["Project Name"])
    # pprint(f"dir(df)={dir(df)}")
    # if len(projects_data) > 0:
    #     df = df.style.sort_values(by='Project Name')
    display(df)
else:
    pprint("Problem in fetching projects")


# In[5]:


project_id = "cb4ef267-0e87-4991-83b6-d5116055db7b" # HPdK_MassPropertiesModelSmall


# # Get branches

# In[6]:


branches_url = f"{host}/projects/{project_id}/branches" 
branches_response = requests.get(branches_url)

if branches_response.status_code == 200:
    branches = branches_response.json()
    branches_data = list(map(lambda b: {'Branch Name':b['name'], 'Branch ID':b['@id'], 'Ref Commit (head)':b['head']}, branches))

    df = pd.DataFrame.from_records(branches_data).sort_values(by='Branch Name').style.hide(axis='index')
    display(df)
else:
    pprint(f"Problem in fetching branches from project {project_id}")
    pprint(branches_response)


# In[7]:


commit_id = "723b31bd-2f6e-4fbe-8fc5-9ebdde0d8d4c"


# # Get PartUsage elements from commit

# In[8]:


element_get_url = f"{host}/projects/{project_id}/commits/{commit_id}/elements" 

element_get_response = requests.get(element_get_url, params={"page[size]":9999})

# if element_get_response.status_code == 200:
if element_get_response.ok:
    elements = element_get_response.json()
    with open("HPdK_MassPropertiesModelSmall_dump.json", "w", encoding="utf-8") as json_file:
        json.dump(elements, json_file)
    part_usages = list(filter(lambda x: x["@type"] == "PartUsage", elements))
    pprint(f"response contains {len(elements)} elements")
    pprint(f"response contains {len(part_usages)} PartUsage elements")
    pprint("PartUsage elements:")
    pprint(sorted(part_usages, key=lambda x: x["declaredShortName"]))
    # elements_data = list(map(lambda b: {'Type':b['@type'], 'Element Name':b['declaredName'], 'Element Short':b['declareShortName'], 'Element ID':b['@id']}, elements))
    # df = pd.DataFrame.from_records(elements_data).sort_values(by=['Type', 'Element Name']).style.hide(axis='index')
    # pprint(elements_data)
    # display(df)
else:
    print(f"response status code={element_get_response.status_code} ({element_get_response.reason})")
    print(f"Problem in fetching elements from project {project_id} at commit {commit_id}.")


# In[ ]:




