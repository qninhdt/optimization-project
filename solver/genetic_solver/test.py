import cppimport

cppimport.settings["release_mode"] = True
cppimport.settings["force_rebuild"] = True

import cppimport.import_hook

import hello

print(hello.square(1))
