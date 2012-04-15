#
# The Sejits enabled version of fftlib.
#

import numpy as np

#
#  I hate this. there has got to be a more portable way to expose
#  the path to the template directory to modules calling this
#  module.
#
import inspect, os
this_file = inspect.currentframe().f_code.co_filename
installDir = os.path.dirname(this_file)

print this_file

class lib(object):
    

    def __init__(self):
        self.pure_python = True

    def fftn_R2C(self, x_in):
        '''Sejits enabled N-dim DFT (complex) of a real array'''

        import asp.codegen.templating.template as template
        tempFile = installDir +  "/templates/fftn_r2c.mako"
        mytemplate = template.Template(filename=tempFile, disable_unicode=True)
        rendered = mytemplate.render(num_dims=3)

        import asp.jit.asp_module as asp_module
        mod = asp_module.ASPModule()
        #
        # Specify function to use from the rendered template
        #
        mod.add_function("SJSfftn_R2C", rendered)
        #
        # Specify libraries for the link phase
        #
        import os
        libs   = ["fftw3","m"]
        incDir = []
        libDir = []
        mod.add_library("FFTW library",incDir,libDir,libs)
        #
        # Add header files needed by the template
        #
        incFile = "fftw3.h"
        mod.add_header(incFile)
        #
        # add numpy libraries used inside the template generated code
        #
	import asp.codegen.cpp_ast as cpp_ast
        mod.add_library("numpy",[np.get_include()+"/numpy"])
        mod.add_header("arrayobject.h")
        #
        # Add inititalization of numpy array package
        #
        mod.add_to_init([cpp_ast.Statement("import_array();")])
        #
        #  Execute the function from within the module and return output
        #
        return mod.SJSfftn_R2C(x_in)

    def ifftn_C2R(self,x_in):
        '''Sejits enabled N-dim inv DFT complex to real'''
        import asp.codegen.templating.template as template
        tempFile = installDir +  "/templates/ifftn_c2r.mako"
        mytemplate = template.Template(filename=tempFile, disable_unicode=True)
        rendered = mytemplate.render(num_dims=3)
        import asp.jit.asp_module as asp_module
        mod = asp_module.ASPModule()
        mod.add_function("SJSifftn_C2R", rendered)
        libDir = []
        libs   = ["fftw3","m"]
        incDir = []
        mod.add_library("FFTW library",incDir,libDir,libs)
        incFile = "fftw3.h"
        mod.add_header(incFile)
	import asp.codegen.cpp_ast as cpp_ast
        mod.add_library("numpy",[np.get_include()+"/numpy"])
        mod.add_header("arrayobject.h")
        mod.add_to_init([cpp_ast.Statement("import_array();")])
        return mod.SJSifftn_C2R(x_in)

