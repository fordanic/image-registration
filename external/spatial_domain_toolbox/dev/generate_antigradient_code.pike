#!/usr/bin/env pike

string reindent(string s)
{
    array(string) result = ({});
    int indent_level = 0;
    array(int) open_parenthesis_positions = ({});
    foreach (s / "\n", string row)
    {
	row = String.trim_whites(row);
	if (sizeof(open_parenthesis_positions) > 0)
	    row = " " * (1 + open_parenthesis_positions[-1]) + row;
	else
	    row = " " * 2 * indent_level + row;

	for (int k = 0; k < sizeof(row); k++)
	{
	    string c = row[k..k];
	    if ((c == "{" || c == "}") && sizeof(open_parenthesis_positions) > 0)
	    {
		werror("syntax error\n");
		exit(1);
	    }
	    else if (c == "{")
		indent_level++;
	    else if (c == "}")
		indent_level--;
	    else if (c == "(")
		open_parenthesis_positions += ({k});
	    else if (c == ")")
		open_parenthesis_positions = open_parenthesis_positions[0..sizeof(open_parenthesis_positions) - 2];
	}

	if (has_value(row, "}"))
	    row = row[2..];
	
	result += ({row});
    }
    
    return result * "\n";
}

string generate_downsample3D()
{
    string s = "";
    
    s += "void\n";
    s += "downsample3D(double *rhs, int M, int N, int P,\n";
    s += "             double *rhs_coarse, int Mhalf, int Nhalf, int Phalf)\n";
    s += "{\n";
    s += "  int i, j, p;\n";
    s += "  int index1;\n";
    s += "  int index2;\n";
    s += "  int MN = M * N;\n";

    for (int oddness = 0; oddness < 8; oddness++)
    {
	int Modd = oddness % 2;
	int Nodd = (oddness / 2) % 2;
	int Podd = (oddness / 4) % 2;
	s += "\n";
	s += sprintf("if (M %% 2 == %d && N %% 2 == %d && P %% 2 == %d)\n",
		     Modd, Nodd, Podd);
	s += "{\n";

	array(string) i_values, j_values, p_values;
	
	if (Modd)
	    i_values = ({"0", "i", "(Mhalf - 1)"});
	else
	    i_values = ({"i"});
	
	if (Nodd)
	    j_values = ({"0", "j", "(Nhalf - 1)"});
	else
	    j_values = ({"j"});
	
	if (Podd)
	    p_values = ({"0", "p", "(Phalf - 1)"});
	else
	    p_values = ({"p"});

	foreach (p_values, string p_value)
	{
	    if (p_value != p_values[0])
		s += "\n";
	    
	    if (p_value == "p")
	    {
		if (Podd)
		    s += "for (p = 1; p < Phalf - 1; p++)\n";
		else
		    s += "for (p = 0; p < Phalf; p++)\n";
		s += "{\n";
	    }

	    foreach (j_values, string j_value)
	    {
		if (j_value != j_values[0])
		    s += "\n";
	    
		if (j_value == "j")
		{
		    if (Nodd)
			s += "for (j = 1; j < Nhalf - 1; j++)\n";
		    else
			s += "for (j = 0; j < Nhalf; j++)\n";
		    s += "{\n";
		}

		foreach (i_values, string i_value)
		{
		    if (i_value != i_values[0])
			s += "\n";
	    
		    if (i_value == "i")
		    {
			if (Modd)
			    s += "for (i = 1; i < Mhalf - 1; i++)\n";
			else
			    s += "for (i = 0; i < Mhalf; i++)\n";
			s += "{\n";
		    }


		    s += sprintf("index1 = (%s * Nhalf + %s) * Mhalf + %s;\n",
				 p_value, j_value, i_value);
		    s += sprintf("index2 = (2 * %s * N + 2 * %s) * M + 2 * %s;\n",
				 p_value, j_value, i_value);

		    
		    mapping(string:float) coeffs = ([]);
		    mapping(string:float) cval = (["-":0.5, "*":1.0, "+":0.5]);
		    foreach ("-*+" / "", string a)
			foreach ("-*+" / "", string b)
			    foreach ("-*+" / "", string c)
				coeffs[a+b+c] = 0.5 * cval[a] * cval[b] * cval[c];

		    if (p_value == "0" || !Podd)
		    {
			foreach (indices(coeffs), string index) {
			    if (index[0..0] == "-")
			    {
				coeffs["+" + index[1..2]] *= 2.0;
				m_delete(coeffs, index);
			    }
			}
		    }
		    else if (p_value != "p")
		    {
			foreach (indices(coeffs), string index) {
			    if (index[0..0] == "+")
			    {
				coeffs["-" + index[1..2]] *= 2.0;
				m_delete(coeffs, index);
			    }
			}
		    }
		    
		    if (j_value == "0" || !Nodd)
		    {
			foreach (indices(coeffs), string index) {
			    if (index[1..1] == "-")
			    {
				coeffs[index[0..0] + "+" + index[2..2]] *= 2.0;
				m_delete(coeffs, index);
			    }
			}
		    }
		    else if (j_value != "j")
		    {
			foreach (indices(coeffs), string index) {
			    if (index[1..1] == "+")
			    {
				coeffs[index[0..0] + "-" + index[2..2]] *= 2.0;
				m_delete(coeffs, index);
			    }
			}
		    }
		    
		    if (i_value == "0" || !Modd)
		    {
			foreach (indices(coeffs), string index) {
			    if (index[2..2] == "-")
			    {
				coeffs[index[0..1] + "+"] *= 2.0;
				m_delete(coeffs, index);
			    }
			}
		    }
		    else if (i_value != "i")
		    {
			foreach (indices(coeffs), string index) {
			    if (index[2..2] == "+")
			    {
				coeffs[index[0..1] + "-"] *= 2.0;
				m_delete(coeffs, index);
			    }
			}
		    }

		    s += "rhs_coarse[index1] = (";
		    mapping(float:array(string)) coeff_terms= ([]);
		    foreach (sort(indices(coeffs)), string index)
		    {
			string index_string = "index2";
			if (index[0..0] == "-")
			    index_string += " - MN";
			else if (index[0..0] == "+")
			    index_string += " + MN";

			if (index[1..1] == "-")
			    index_string += " - M";
			else if (index[1..1] == "+")
			    index_string += " + M";
			
			if (index[2..2] == "-")
			    index_string += " - 1";
			else if (index[2..2] == "+")
			    index_string += " + 1";

			if (!coeff_terms[coeffs[index]])
			    coeff_terms[coeffs[index]] = ({index_string});
			else
			    coeff_terms[coeffs[index]] += ({index_string});
		    }
		    
		    int first_coeff = 1;
		    foreach (reverse(sort(indices(coeff_terms))),
			     float coeff_value)
		    {
			string terms = sprintf("rhs[%s]", coeff_terms[coeff_value][*]) * "\n+ ";
			if (coeff_value != 1.0)
			{
			    if (first_coeff)
				first_coeff = 0;
			    else
				s += "+ ";
	    
			    s += sprintf("%s * (%s)\n", format_coeff(coeff_value), terms);
			}
			else
			    s += sprintf("%s\n", terms);
		    }

		    s = s[..sizeof(s)-2];
		    s += ");\n";
		    
		    if (i_value == "i")
			s += "}\n";
		}	    
		
		if (j_value == "j")
		    s += "}\n";
	    }	    
		
	    if (p_value == "p")
		s += "}\n";
	}
	
	s += "}\n";
    }
    s += "}\n";
    return reindent(s);
}

mapping(string:float) combine_coeffs(mapping(string:float) coeffs,
				     array(float) dist)
{
    mapping(string:float) new_coeffs = ([]);
    foreach (indices(coeffs), string pos)
    {
	if (dist[0] > 0.0)
	    new_coeffs[pos + "-"] = coeffs[pos] * dist[0];
	if (dist[1] > 0.0)
	    new_coeffs[pos + "*"] = coeffs[pos] * dist[1];
	if (dist[2] > 0.0)
	    new_coeffs[pos + "+"] = coeffs[pos] * dist[2];
    }

    return new_coeffs;
}

mapping(string:float) mirror_coeffs(mapping(string:float) coeffs, int pos,
				    string from)
{
    string to = "*";
    mapping(string:float) new_coeffs = ([]);
    foreach (coeffs; string index; float value)
    {
	if (index[pos..pos] == from)
	    index = index[..pos-1] + to + index[pos+1..];
	if (!new_coeffs[index])
	    new_coeffs[index] = value;
	else
	    new_coeffs[index] += value;
    }

    return new_coeffs;
}

string format_coeff(float x)
{
    string s = sprintf("%8.6f", x);
    while (s[sizeof(s)-1] == '0')
	s = s[..sizeof(s)-2];
    return s;
}

string do_interpolate(mapping(string:float) coeffs, string out_index)
{
    string s = "";
    s += sprintf("f_out[%s] += (", out_index);

    mapping(float:array(string)) coeff_terms= ([]);
    foreach (sort(indices(coeffs)), string index)
    {
	string index_string = "index1";
	if (index[0..0] == "-")
	    index_string += " - MNhalf";
	else if (index[0..0] == "+")
	    index_string += " + MNhalf";
	
	if (index[1..1] == "-")
	    index_string += " - Mhalf";
	else if (index[1..1] == "+")
	    index_string += " + Mhalf";
	
	if (index[2..2] == "-")
	    index_string += " - 1";
	else if (index[2..2] == "+")
	    index_string += " + 1";
	
	if (!coeff_terms[coeffs[index]])
	    coeff_terms[coeffs[index]] = ({index_string});
	else
	    coeff_terms[coeffs[index]] += ({index_string});
    }

    int first_coeff = 1;
    foreach (reverse(sort(indices(coeff_terms))),
	     float coeff_value)
    {
	string terms = sprintf("v[%s]", coeff_terms[coeff_value][*]) * "\n+ ";
	if (coeff_value != 1.0)
	{
	    if (first_coeff)
		first_coeff = 0;
	    else
		s += "+ ";
	    
	    if (sizeof(coeff_terms[coeff_value]) > 1)
		s += sprintf("%s * (%s)\n", format_coeff(coeff_value), terms);
	    else
		s += sprintf("%s * %s\n", format_coeff(coeff_value), terms);
	}
	else
	    s += sprintf("%s\n", terms);
    }
    
    s = s[..sizeof(s)-2];
    s += ");\n";
    return s;
}

string split_in_cases(mapping(string:float) coeffs, int pos,
		      array(string) vars, array(string) upper_limits,
		      array(int) skip_upper_limits, string out_index)
{
    string s = "";
    
    if (has_value(column(indices(coeffs), pos), '-'))
    {
	s += sprintf("if (%s == 0)\n", vars[pos]);
	s += "{\n";
	mapping(string:float) new_coeffs = mirror_coeffs(coeffs, pos, "-"); 
	if (pos < sizeof(vars) - 1)
	    s += split_in_cases(new_coeffs, pos + 1, vars, upper_limits,
				skip_upper_limits, out_index);
	else
	    s += do_interpolate(new_coeffs, out_index);
	s += "}\n";
	s += "else\n";
	s += "{\n";
	if (pos < sizeof(vars) - 1)
	    s += split_in_cases(coeffs, pos + 1, vars, upper_limits,
				skip_upper_limits, out_index);
	else
	    s += do_interpolate(coeffs, out_index);
	s += "}\n";
    }
    else if (has_value(column(indices(coeffs), pos), '+'))
    {
	s += sprintf("if (%s < %s - 1)\n", vars[pos], upper_limits[pos]);
	s += "{\n";
	if (pos < sizeof(vars) - 1)
	    s += split_in_cases(coeffs, pos + 1, vars, upper_limits,
				skip_upper_limits, out_index);
	else
	    s += do_interpolate(coeffs, out_index);
	s += "}\n";
	if (!skip_upper_limits[pos])
	{
	    s += "else\n";
	    s += "{\n";
	    mapping(string:float) new_coeffs = mirror_coeffs(coeffs, pos, "+"); 
	    if (pos < sizeof(vars) - 1)
		s += split_in_cases(new_coeffs, pos + 1, vars, upper_limits,
				    skip_upper_limits, out_index);
	    else
		s += do_interpolate(new_coeffs, out_index);
	    s += "}\n";
	}
    }
    else
    {
	if (pos < sizeof(vars) - 1)
	    s += split_in_cases(coeffs, pos + 1, vars, upper_limits,
				skip_upper_limits, out_index);
	else
	    s += do_interpolate(coeffs, out_index);
    }

    return s;
}

string generate_upsample3D()
{
    string s = "";
    
    s += "void\n";
    s += "upsample3D(double *rhs, int M, int N, int P,\n";
    s += "           double *v, int Mhalf, int Nhalf, int Phalf,\n";
    s += "           double *f_out)\n";
    s += "{\n";
    s += "  int i, j, p;\n";
    s += "  int index1;\n";
    s += "  int index2;\n";
    s += "  int MN = M * N;\n";
    s += "  int MNhalf = Mhalf * Nhalf;\n";

    for (int oddness = 0; oddness < 8; oddness++)
    {
	int Modd = oddness % 2;
	int Nodd = (oddness / 2) % 2;
	int Podd = (oddness / 4) % 2;
	s += "\n";
	s += sprintf("if (M %% 2 == %d && N %% 2 == %d && P %% 2 == %d)\n",
		     Modd, Nodd, Podd);
	s += "{\n";

	array(string) i_values, j_values, p_values;
	
	if (Modd)
	    i_values = ({"0", "i", "(Mhalf - 1)"});
	else
	    i_values = ({"i"});
	
	if (Nodd)
	    j_values = ({"0", "j", "(Nhalf - 1)"});
	else
	    j_values = ({"j"});
	
	if (Podd)
	    p_values = ({"0", "p", "(Phalf - 1)"});
	else
	    p_values = ({"p"});


	s += "for (p = 0; p < Phalf; p++)\n";
	s += "{\n";
	s += "for (j = 0; j < Nhalf; j++)\n";
	s += "{\n";
	s += "for (i = 0; i < Mhalf; i++)\n";
	s += "{\n";

	s += "index1 = (p * Nhalf + j) * Mhalf + i;\n";
	s += "index2 = (2 * p * N + 2 * j) * M + 2 * i;\n";

	for (int k = 0; k < 8; k++)
	{
	    int i = k % 2;
	    int j = (k / 2) % 2;
	    int p = k / 4;
	    mapping (string:float) coeffs = (["" : 1.0]);

	    if (k > 0)
		s += "\n";
	    
	    if (Podd && p == 0)
		coeffs = combine_coeffs(coeffs, ({0.0, 1.0, 0.0}));
	    else if (Podd && p == 1)
		coeffs = combine_coeffs(coeffs, ({0.0, 0.5, 0.5}));
	    else if (!Podd && p == 0)
		coeffs = combine_coeffs(coeffs, ({0.25, 0.75, 0.0}));
	    else if (!Podd && p == 1)
		coeffs = combine_coeffs(coeffs, ({0.0, 0.75, 0.25}));
	    
	    if (Nodd && j == 0)
		coeffs = combine_coeffs(coeffs, ({0.0, 1.0, 0.0}));
	    else if (Nodd && j == 1)
		coeffs = combine_coeffs(coeffs, ({0.0, 0.5, 0.5}));
	    else if (!Nodd && j == 0)
		coeffs = combine_coeffs(coeffs, ({0.25, 0.75, 0.0}));
	    else if (!Nodd && j == 1)
		coeffs = combine_coeffs(coeffs, ({0.0, 0.75, 0.25}));
	    
	    if (Modd && i == 0)
		coeffs = combine_coeffs(coeffs, ({0.0, 1.0, 0.0}));
	    else if (Modd && i == 1)
		coeffs = combine_coeffs(coeffs, ({0.0, 0.5, 0.5}));
	    else if (!Modd && i == 0)
		coeffs = combine_coeffs(coeffs, ({0.25, 0.75, 0.0}));
	    else if (!Modd && i == 1)
		coeffs = combine_coeffs(coeffs, ({0.0, 0.75, 0.25}));

	    string out_index = "index2";
	    if (p == 1)
		out_index += " + MN";
	    if (j == 1)
		out_index += " + M";
	    if (i == 1)
		out_index += " + 1";
	    
	    s += split_in_cases(coeffs, 0, ({"p", "j", "i"}),
				({"Phalf", "Nhalf", "Mhalf"}),
				({p && Podd, j && Nodd, i && Modd}),
				out_index);
	}
	
	s += "}\n";
	s += "}\n";
	s += "}\n";
	s += "}\n";
    }
    
    s += "}\n";
    
    return reindent(s);
}

int main(int argc, array(string) argv)
{
    write(generate_downsample3D());
    write("\n\n");
    write(generate_upsample3D());
}
