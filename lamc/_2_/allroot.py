import os
for root, dirs, files in os.walk("."):  
    for filename in files:
        if ".hist" in filename:
            os.system(f"h2root {filename}")
            os.system(f".q")

