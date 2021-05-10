import json  
       
# Data to be written  
dictionary ={  
  "id": "04",  
  "name": "sunil",  
  "depatment": "HR"
}  
       
# Serializing json   
json_object = json.dumps(dictionary, indent = 2)  
print(json_object) 