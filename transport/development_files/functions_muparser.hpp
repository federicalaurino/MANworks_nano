#ifndef M3D1D_MUPARSER_FUNCTIONS_HPP_
#define M3D1D_MUPARSER_FUNCTIONS_HPP_

 namespace getfem {
double CC1=0.0, kk1=0.0, RR1=0.0;
std::string expr1="";


double generic_function(const bgeot::base_node & x){

	//Define p
	mu::Parser p;
	//Define variables
		/* We define both coordinates x,y,z,
		 * which are needed for the function,
		 * and other useful parameters, namely C,R,k,
		 * that you can pass in the .param file.
		 * Note that you can add in this file any other 
		 * variable you want. */ 
	double xx=x[0], yy=x[1], zz=x[2];
	p.DefineVar("x", &xx);
	p.DefineVar("y", &yy);
	p.DefineVar("z", &zz);
	p.DefineVar("C", &CC1);	
	p.DefineVar("k", &kk1);
	p.DefineVar("R", &RR1);
	//p.DefineVar("new variable", double variable);	
	
	//Define the expression
		/* The string should be already modified
		 * in the code */
	p.SetExpr(expr1);
	return p.Eval();
}





} // end namespace getfem
#endif
