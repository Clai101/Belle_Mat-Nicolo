import os
os.system("mkdir root")
os.system("mkdir hist")
for root, dirs, files in os.walk("."):  
    for filename in files:
        if ".hist" in filename:
            os.system(f"h2root {filename}")
os.system(r"mv *.hist hist")
os.system(r"mv *.root root")

