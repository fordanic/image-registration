#!/usr/bin/env pike

string reindent(string s)
{
  array(string) result = ({});
  int indent_level = 0;
  array(int) open_parenthesis_positions = ({});
  int previous_line_if_statement = 0;
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

    if (!has_value(row, "{") && previous_line_if_statement)
      row = "  " + row;
    
    result += ({row});

    previous_line_if_statement = (has_value(row, "if (")
				  || has_value(row, "else"));
  }
  
  return result * "\n";
}

string generate_downsample3D()
{
  string s = "";
  
  s += "static void\n";
  s += "downsample3D(double *rhs, int M, int N, int P,\n";
  s += "             double *rhs_coarse, int Mhalf, int Nhalf, int Phalf,\n";
  s += "             double *weight, double *coarse_weight)\n";
  s += "{\n";
  s += "  int i, j, p;\n";
  s += "  int index1;\n";
  s += "  int index2;\n";
  s += "  int MN = M * N;\n";
  s += "  double w[3][3][3];\n";
  s += "  double sum;\n";
  
  for (int oddness = 0; oddness < 8; oddness++)
  {
    int Modd = oddness % 2;
    int Nodd = (oddness / 2) % 2;
    int Podd = (oddness / 4) % 2;
    s += "\n";
    s += sprintf("if (M %% 2 == %d && N %% 2 == %d && P %% 2 == %d)\n",
		 Modd, Nodd, Podd);
    s += "{\n";
    
    s += "    for (p = 0; p < Phalf; p++)\n";
    s += "    {\n";
    s += "      for (j = 0; j < Nhalf; j++)\n";
    s += "      {\n";
    s += "        for (i = 0; i < Mhalf; i++)\n";
    s += "        {\n";
    s += "          index1 = (p * Nhalf + j) * Mhalf + i;\n";
    s += "          index2 = (2 * p * N + 2 * j) * M + 2 * i;\n";
    s += "\n";
    
    array(string) i_conds, j_conds, p_conds;
    array(string) i_offsets, j_offsets, p_offsets;
    array(float) i_coeffs, j_coeffs, p_coeffs;
    
    if (Modd)
    {
      i_conds  = ({"i > 0", "", "i < Mhalf - 1"});
      i_coeffs = ({0.5, 1.0, 0.5});
    }
    else
    {
      i_conds  = ({"", "", ""});
      i_coeffs = ({0.0, 1.0, 1.0});
    }
    i_offsets = ({" - 1", "", " + 1"});
    
    if (Nodd)
    {
      j_conds  = ({"j > 0", "", "j < Nhalf - 1"});
      j_coeffs = ({0.5, 1.0, 0.5});
    }
    else
    {
      j_conds  = ({"", "", ""});
      j_coeffs = ({0.0, 1.0, 1.0});
    }
    j_offsets = ({" - M", "", " + M"});
    
    if (Podd)
    {
      p_conds  = ({"p > 0", "", "p < Phalf - 1"});
      p_coeffs = ({0.5, 1.0, 0.5});
    }
    else
    {
      p_conds  = ({"", "", ""});
      p_coeffs = ({0.0, 1.0, 1.0});
    }
    p_offsets = ({" - MN", "", " + MN"});
    
    string s2 = "";
    array(string) terms = ({});
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
	for (int p = 0; p < 3; p++)
	{
	  float coeff = i_coeffs[i] * j_coeffs[j] * p_coeffs[p];
	  if (coeff == 0.0)
	    continue;
	  
	  array(string) conds = ({i_conds[i], j_conds[j],
				  p_conds[p]}) - ({""});
	  
	  string w = sprintf("w[%d][%d][%d]", i, j, p);
	    
	  string index = sprintf("index2%s%s%s", p_offsets[p],
				 j_offsets[j], i_offsets[i]);
	  
	  if (coeff == 1.0)
	    s += sprintf("%s = weight[%s];\n", w, index);
	  else
	    s += sprintf("%s = %s * VAL(%s, weight[%s]);\n",
			 w, format_coeff(coeff), conds * " && ", index);
	  
	  s2 += sprintf("if (%s > 0)\n", w);
	  s2 += sprintf("result += %s * rhs[%s];\n", w, index);
	  
	  terms += ({w});
	}

    s += "\n";
    s += sprintf("sum = %s;\n", terms * " + ");
    s += "coarse_weight[index1] = sum;\n";
    s += "\n";
    s += "if (sum > 0)\n";
    s += "{\n";
    s += "double result = 0;\n";
    s += s2;
    s += "\n";
    s += "rhs_coarse[index1] = 8 / sum * result;\n";
    s += "}\n";
    s += "}\n";
    s += "}\n";
    s += "}\n";
    s += "}\n";
  }
  s += "}\n";
  return reindent(s);
}


string generate_galerkin3D()
{
  string s = "";
  
  s += "static void\n";
  s += "galerkin3D(int level, int M, int N, int P, int Mhalf, int Nhalf, int Phalf,\n";
  s += "           double *weight, double *coarse_weight)\n";
  s += "{\n";
  s += "  int i, j, p;\n";
  s += "  int k;\n";
  s += "  double *lhs;\n";
  s += "  double *lhs_coarse;\n";
  s += "  int MN = M * N;\n";
  s += "  int MNhalf = Mhalf * Nhalf;\n";
  s += "  double lw[3][3][3];\n";
  s += "  double mean;\n";
  s += "\n";
  s += "  if (data.lhs[level + 1] != NULL)\n";
  s += "    return;\n";
  s += "\n";
  s += "  data.lhs[level + 1] = mxCalloc(27 * Mhalf * Nhalf * Phalf,\n";
  s += "				 sizeof(*data.lhs[level + 1]));\n";
  s += "  lhs = data.lhs[level];\n";
  s += "  lhs_coarse = data.lhs[level + 1];\n";
  s += "\n";
  s += "  for (p = 0; p < Phalf; p++)\n";
  s += "  {\n";
  s += "    for (j = 0; j < Nhalf; j++)\n";
  s += "    {\n";
  s += "      for (i = 0; i < Mhalf; i++)\n";
  s += "      {\n";
  s += "        int index1 = ((p * Nhalf + j) * Mhalf + i);\n";
  s += "        int index2 = ((2 * p * N + 2 * j) * M + 2 * i);\n";
  s += "        double stencil1[3][3][3];\n";
  s += "        double stencil2[5][5][5];\n";
  s += "        double stencil3[3][3][3];\n";
  s += "        double mask1[3][3][3];\n";
  s += "        int u, v, w;\n";
  s += "\n";
  s += "        for (u = 0; u < 5; u++)\n";
  s += "        {\n";
  s += "          for (v = 0; v < 5; v++)\n";
  s += "          {\n";
  s += "            for (w = 0; w < 5; w++)\n";
  s += "            {\n";
  s += "              stencil2[u][v][w] = 0;\n";
  s += "              if (u < 3 && v < 3 && w < 3)\n";
  s += "              {\n";
  s += "                stencil1[u][v][w] = 0;\n";
  s += "                stencil3[u][v][w] = 0;\n";
  s += "              }\n";
  s += "            }\n";
  s += "          }\n";
  s += "        }\n";
  
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int p = 0; p < 3; p++)
      {
	string index = sprintf("index1%s%s%s",
			       ({" - MNhalf", "", " + MNhalf"})[p],
			       ({" - Mhalf", "", " + Mhalf"})[j],
			       ({" - 1", "", " + 1"})[i]);

	array(string) conditions = ({});
	conditions += ({"i > 0", "", "i < Mhalf - 1"})[i..i];
	conditions += ({"j > 0", "", "j < Nhalf - 1"})[j..j];
	conditions += ({"p > 0", "", "p < Phalf - 1"})[p..p];
	conditions -= ({""});

	s += sprintf("mask1[%d][%d][%d] = ", i, j, p);
	
	if (sizeof(conditions) == 0)
	  s += sprintf("coarse_weight[%s];\n", index);
	else
	  s += sprintf("VAL(%s,\n coarse_weight[%s]);\n",
		       conditions * " && ", index);
      }
  
  for (int oddness = 0; oddness < 8; oddness++)
  {
    int Modd = oddness % 2;
    int Nodd = (oddness / 2) % 2;
    int Podd = (oddness / 4) % 2;
    s += "\n";
    s += sprintf("if (M %% 2 == %d && N %% 2 == %d && P %% 2 == %d)\n",
		 Modd, Nodd, Podd);
    s += "{\n";
    
    array(string) i_conds, j_conds, p_conds;
    array(string) i_offsets, j_offsets, p_offsets;
    array(float) i_coeffs, j_coeffs, p_coeffs;
    
    if (Modd)
    {
      i_conds  = ({"i > 0", "", "i < Mhalf - 1"});
      i_coeffs = ({0.5, 1.0, 0.5});
    }
    else
    {
      i_conds  = ({"", "", ""});
      i_coeffs = ({0.0, 1.0, 1.0});
    }
    i_offsets = ({" - 1", "", " + 1"});
    
    if (Nodd)
    {
      j_conds  = ({"j > 0", "", "j < Nhalf - 1"});
      j_coeffs = ({0.5, 1.0, 0.5});
    }
    else
    {
      j_conds  = ({"", "", ""});
      j_coeffs = ({0.0, 1.0, 1.0});
    }
    j_offsets = ({" - M", "", " + M"});
    
    if (Podd)
    {
      p_conds  = ({"p > 0", "", "p < Phalf - 1"});
      p_coeffs = ({0.5, 1.0, 0.5});
    }
    else
    {
      p_conds  = ({"", "", ""});
      p_coeffs = ({0.0, 1.0, 1.0});
    }
    p_offsets = ({" - MN", "", " + MN"});
    
    string s2 = "";
    array(string) terms = ({});
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
	for (int p = 0; p < 3; p++)
	{
	  float coeff = i_coeffs[i] * j_coeffs[j] * p_coeffs[p];
	  if (coeff == 0.0)
	    continue;
	  
	  array(string) conds = ({i_conds[i], j_conds[j],
				  p_conds[p]}) - ({""});
	  
	  string w = sprintf("lw[%d][%d][%d]", i, j, p);
	    
	  string index = sprintf("index2%s%s%s", p_offsets[p],
				 j_offsets[j], i_offsets[i]);
	  
	  if (coeff == 1.0)
	    s += sprintf("%s = weight[%s];\n", w, index);
	  else
	    s += sprintf("%s = %s * VAL(%s, weight[%s]);\n",
			 w, format_coeff(coeff), conds * " && ", index);
	  
	  s2 += sprintf("stencil1[%d][%d][%d] = %s / mean;\n", i, j, p, w);
	  
	  terms += ({w});
	}

    s += "\n";
    s += sprintf("mean = (%s) / 8;\n", terms * " + ");
    s += "\n";
    s += "if (mean == 0)\n";
    s += "continue;\n";
    s += "\n";
    s += s2;
    s += "}\n";
  }

  s += "\n";
  s += "for (u = 0; u < 3; u++)\n";
  s += "{\n";
  s += "for (v = 0; v < 3; v++)\n";
  s += "{\n";
  s += "for (w = 0; w < 3; w++)\n";
  s += "{\n";
  s += "if (stencil1[u][v][w] != 0)\n";
  s += "{\n";
  s += "int index = 27 * (index2 + ((w - 1) * N + v - 1) * M + u - 1);\n";
  s += "if (lhs[index + 13] != 0)\n";
  s += "{\n";
  s += "int k;\n";
  s += "for (k = 0; k < 27; k++)\n";
  s += "{\n";
  s += "  int a = (k % 3);\n";
  s += "  int b = ((k / 3) % 3);\n";
  s += "  int c = ((k / 9) % 3);\n";
  s += "  stencil2[u + a][v + b][w + c] += stencil1[u][v][w] * lhs[index + k];\n";
  s += "}\n";
  s += "}\n";
  s += "}\n";
  s += "}\n";
  s += "}\n";
  s += "}\n";

  for (int oddness = 0; oddness < 8; oddness++)
  {
    int Modd = oddness % 2;
    int Nodd = (oddness / 2) % 2;
    int Podd = (oddness / 4) % 2;
    s += "\n";
    s += sprintf("if (M %% 2 == %d && N %% 2 == %d && P %% 2 == %d)\n",
		 Modd, Nodd, Podd);
    s += "{\n";
    
    string ihalf, jhalf, phalf;
    
    s += sprintf("for (u = %d; u < 5; u++)\n", 1 - Modd);
    s += "{\n";
    s += sprintf("for (v = %d; v < 5; v++)\n", 1 - Nodd);
    s += "{\n";
    s += sprintf("for (w = %d; w < 5; w++)\n", 1 - Podd);
    s += "{\n";
    s += "double alpha1, alpha2, alpha3;\n";
    s += "double unw, dnw, une, dne, usw, dsw, use, dse;\n";
    s += "int uu, vv, ww;\n";
    s += "double sum;\n";
    s += "\n";
    s += "if (u % 2 == 0)\n";
    s += sprintf("alpha1 = %s;\n", Modd ? "0" : "0.75");
    s += "else\n";
    s += sprintf("alpha1 = %s;\n", Modd ? "0.5" : "0.25");
    s += "\n";
    s += "if (v % 2 == 0)\n";
    s += sprintf("alpha2 = %s;\n", Nodd ? "0" : "0.75");
    s += "else\n";
    s += sprintf("alpha2 = %s;\n", Nodd ? "0.5" : "0.25");
    s += "\n";
    s += "if (w % 2 == 0)\n";
    s += sprintf("alpha3 = %s;\n", Podd ? "0" : "0.75");
    s += "else\n";
    s += sprintf("alpha3 = %s;\n", Podd ? "0.5" : "0.25");
    s += "\n";

    s += sprintf("uu = %s / 2;\n", Modd ? "u" : "(u - 1)");
    s += sprintf("vv = %s / 2;\n", Nodd ? "v" : "(v - 1)");
    s += sprintf("ww = %s / 2;\n", Podd ? "w" : "(w - 1)");
    
    s += sprintf("%s;\n", upsample_w2(0, 0, 0));
    s += sprintf("%s;\n", upsample_w2(0, 0, 1));
    s += sprintf("%s;\n", upsample_w2(0, 1, 0));
    s += sprintf("%s;\n", upsample_w2(0, 1, 1));
    s += sprintf("%s;\n", upsample_w2(1, 0, 0));
    s += sprintf("%s;\n", upsample_w2(1, 0, 1));
    s += sprintf("%s;\n", upsample_w2(1, 1, 0));
    s += sprintf("%s;\n", upsample_w2(1, 1, 1));
    s += "\n";
    
    s += "sum = unw + dnw + une + dne + usw + dsw + use + dse;\n";
    s += "\n";
    
    s += "if (sum > 0)\n";
    s += "{\n";
    s += upsample_contribution2(0, 0, 0);
    s += upsample_contribution2(0, 0, 1);
    s += upsample_contribution2(0, 1, 0);
    s += upsample_contribution2(0, 1, 1);
    s += upsample_contribution2(1, 0, 0);
    s += upsample_contribution2(1, 0, 1);
    s += upsample_contribution2(1, 1, 0);
    s += upsample_contribution2(1, 1, 1);
    s += "}\n";
    s += "}\n";
    s += "}\n";
    s += "}\n";
    s += "}\n";
  }

  s += "\n";
  s += "for (k = 0; k < 27; k++)\n";
  s += "{\n";
  s += "  int a = (k % 3);\n";
  s += "  int b = ((k / 3) % 3);\n";
  s += "  int c = ((k / 9) % 3);\n";
  s += "  lhs_coarse[27 * index1 + k] = stencil3[a][b][c];\n";
  s += "}\n";
  s += "\n";

  s += "for (u = 0; u < 27; u++)\n";
  s += "{\n";
  s += "  if (u != 13 && lhs_coarse[27 * index1 + u] < 0)\n";
  s += "  {\n";
  s += "    lhs_coarse[27 * index1 + 13] += lhs_coarse[27 * index1 + u];\n";
  s += "    lhs_coarse[27 * index1 + u] = 0;\n";
  s += "  }\n";
  s += "}\n";
  s += "}\n";
  s += "}\n";
  s += "}\n";
  s += "}\n";
  return reindent(s);
}


string format_coeff(float x)
{
    string s = sprintf("%8.6f", x);
    while (s[sizeof(s)-1] == '0')
	s = s[..sizeof(s)-2];
    return s;
}


string upsample_w(int i, int j, int p, int Modd, int Nodd, int Podd)
{
  string s = sprintf("%s%s%s = (", p ? "d" : "u", i ? "s" : "n", j ? "e" : "w");
  s += i ? "alpha1" : "(1 - alpha1)";
  s += " * ";
  s += j ? "alpha2" : "(1 - alpha2)";
  s += " * ";
  s += p ? "alpha3" : "(1 - alpha3)";
  s += " * \n";

  string index = sprintf("index2%s%s%s", p ? " + MNhalf" : "",
			 j ? " + Mhalf" : "", i ? " + 1" : "");

  array(string) conditions = ({});
  if (i)
    conditions += ({"i < Mhalf - 1"});
  else if (!Modd)
    conditions += ({"i > 0"});
  
  if (j)
    conditions += ({"j < Nhalf - 1"});
  else if (!Nodd)
    conditions += ({"j > 0"});
  
  if (p)
    conditions += ({"p < Phalf - 1"});
  else if (!Podd)
    conditions += ({"p > 0"});
  
  if (sizeof(conditions) == 0)
    s += sprintf("coarse_weight[%s])", index);
  else
    s += sprintf("VAL(%s,\n coarse_weight[%s]))", conditions * " && ", index);

  return s;
}

string upsample_w2(int i, int j, int p)
{
  string s = sprintf("%s%s%s = (", p ? "d" : "u", i ? "s" : "n", j ? "e" : "w");
  s += i ? "alpha1" : "(1 - alpha1)";
  s += " * ";
  s += j ? "alpha2" : "(1 - alpha2)";
  s += " * ";
  s += p ? "alpha3" : "(1 - alpha3)";
  s += " * \n";

  array(string) conditions = ({});
  if (i)
    conditions += ({"alpha1 > 0"});
  if (j)
    conditions += ({"alpha2 > 0"});
  if (p)
    conditions += ({"alpha3 > 0"});

  string mask1 = sprintf("mask1[uu%s][vv%s][ww%s]",
			 i ? " + 1" : "",
			 j ? " + 1" : "",
			 p ? " + 1" : "");
  
  if (sizeof(conditions) == 0)
    s += mask1 + ")";
  else
    s += sprintf("VAL(%s,\n%s))", conditions * " && ", mask1);

  return s;
}

string upsample_contribution(int i, int j, int p)
{
  string w = sprintf("%s%s%s", p ? "d" : "u", i ? "s" : "n", j ? "e" : "w");
  string index = sprintf("index2%s%s%s", p ? " + MNhalf" : "",
			 j ? " + Mhalf" : "", i ? " + 1" : "");
  string s = sprintf("if (%s > 0)\n", w);
  s += sprintf("contribution += %s * v[%s];\n", w, index);
  return s;
}

string upsample_contribution2(int i, int j, int p)
{
  string w = sprintf("%s%s%s", p ? "d" : "u", i ? "s" : "n", j ? "e" : "w");
  string index = sprintf("[uu%s][vv%s][ww%s]",
			 i ? " + 1" : "",
			 j ? " + 1" : "",
			 p ? " + 1" : "");
  string s = sprintf("if (%s > 0)\n", w);
  s += sprintf("stencil3%s += %s * stencil2[u][v][w] / sum;\n", index, w);
  return s;
}

string generate_upsample3D()
{
  string s = "";
    
  s += "static void\n";
  s += "upsample3D(double *rhs, int M, int N, int P,\n";
  s += "           double *v, int Mhalf, int Nhalf, int Phalf,\n";
  s += "           double *f_out, double *weight, double *coarse_weight)\n";
  s += "{\n";
  s += "  int i, j, p;\n";
  s += "  int index1;\n";
  s += "  int index2;\n";
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
    
    string ihalf, jhalf, phalf;
    
    if (Modd)
      ihalf = "i / 2";
    else
      ihalf = "(i - 1) / 2";
    
    if (Nodd)
      jhalf = "j / 2";
    else
      jhalf = "(j - 1) / 2";
    
    if (Podd)
      phalf = "p / 2";
    else
      phalf = "(p - 1) / 2";
    
    s += "for (p = 0; p < P; p++)\n";
    s += "{\n";
    s += "for (j = 0; j < N; j++)\n";
    s += "{\n";
    s += "for (i = 0; i < M; i++)\n";
    s += "{\n";
    s += "double alpha1, alpha2, alpha3;\n";
    s += "double unw, dnw, une, dne, usw, dsw, use, dse;\n";
    s += "double sum;\n";
    s += "\n";
    s += "index1 = (p * N + j) * M + i;\n";
    s += sprintf("index2 = ((%s) * Nhalf + %s) * Mhalf + %s;\n",
		 phalf, jhalf, ihalf);
    s += "\n";
    s += "if (weight[index1] == 0)\n";
    s += "continue;\n";
    s += "\n";
    s += "if (i % 2 == 0)\n";
    s += sprintf("alpha1 = %s;\n", Modd ? "0" : "0.75");
    s += "else\n";
    s += sprintf("alpha1 = %s;\n", Modd ? "0.5" : "0.25");
    s += "\n";
    s += "if (j % 2 == 0)\n";
    s += sprintf("alpha2 = %s;\n", Nodd ? "0" : "0.75");
    s += "else\n";
    s += sprintf("alpha2 = %s;\n", Nodd ? "0.5" : "0.25");
    s += "\n";
    s += "if (p % 2 == 0)\n";
    s += sprintf("alpha3 = %s;\n", Podd ? "0" : "0.75");
    s += "else\n";
    s += sprintf("alpha3 = %s;\n", Podd ? "0.5" : "0.25");
    s += "\n";
    
    s += sprintf("%s;\n", upsample_w(0, 0, 0, Modd, Nodd, Podd));
    s += sprintf("%s;\n", upsample_w(0, 0, 1, Modd, Nodd, Podd));
    s += sprintf("%s;\n", upsample_w(0, 1, 0, Modd, Nodd, Podd));
    s += sprintf("%s;\n", upsample_w(0, 1, 1, Modd, Nodd, Podd));
    s += sprintf("%s;\n", upsample_w(1, 0, 0, Modd, Nodd, Podd));
    s += sprintf("%s;\n", upsample_w(1, 0, 1, Modd, Nodd, Podd));
    s += sprintf("%s;\n", upsample_w(1, 1, 0, Modd, Nodd, Podd));
    s += sprintf("%s;\n", upsample_w(1, 1, 1, Modd, Nodd, Podd));
    s += "\n";
    
    s += "sum = unw + dnw + une + dne + usw + dsw + use + dse;\n";
    s += "\n";
    
    s += "if (sum > 0)\n";
    s += "{\n";
    s += "double contribution = 0;\n";
    s += "\n";
    s += upsample_contribution(0, 0, 0);
    s += upsample_contribution(0, 0, 1);
    s += upsample_contribution(0, 1, 0);
    s += upsample_contribution(0, 1, 1);
    s += upsample_contribution(1, 0, 0);
    s += upsample_contribution(1, 0, 1);
    s += upsample_contribution(1, 1, 0);
    s += upsample_contribution(1, 1, 1);
    s += "\n";
    s += "f_out[index1] += contribution / sum;\n";
    s += "}\n";
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
    write(generate_galerkin3D());
    write("\n\n");
    write(generate_upsample3D());
}
