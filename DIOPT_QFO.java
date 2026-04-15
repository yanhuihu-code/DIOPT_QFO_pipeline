// File: DIOPT pipeline
// Author:	Claire Hu

//==================================================================//
// This program is responsible for process ortholog prediction results
// from various algorithms into an uniform format
//==================================================================//

import java.io.*;
import java.net.*;
import java.lang.*;
import java.util.*;
import java.sql.*;

public class DIOPT_QFO
{
///// step1: assemble gene information based on NCBI gene file (https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz), need to upzip first before run the program
 	public static void assembleGeneInfo(String input, String taxids, String output)
    {
		String strr="";
		String str2 = "";
        try
        {
			File ta = new File(taxids);
			String[] tax1 = new String[100];
			String[] tax2 = new String[100];
			String[] db = new String[100];
			int lt = fillStringArray(ta,tax1,tax2,db);
			System.out.println("species selected:"+lt);
            FileWriter out = new FileWriter(output);
            out.write("geneid\tspeciesid\tsymbol\tNdescription\tlocus_tag\tspecies_specific_geneid\tspecies_specific_geneid_type\tchromosome\tmap_location\tgene_type\n");
            File f = new File(input);
        	BufferedReader b1 = new BufferedReader(new FileReader(f));
        	int fcount = 0;
        	int count_all = 0;
        	int count = 0;
        	String geneid ="";
            String taxid = "";
			String symbol = "";
			String locustag = "";
			String nick = "";
			String dbxrefs = "";
			String chr = "";
			String map = "";
			String desc = "";
			String type = "";
			String xdb = "";
			String xdb_final = "";
			String xdbid = "";
			String tax_final="";
			String tmp1 = "";
			String tmp2 = "";
			boolean want = false;
         	while( (strr=b1.readLine())!=null)//read file line by line
        	{
				if (strr.length()==0 || strr.indexOf("NEWENTRY")>-1) continue;
				count_all++;
          		StringTokenizer st = new StringTokenizer(strr,"\t");
                if (st.countTokens()>=10)
                {
					taxid = st.nextToken();
					geneid = st.nextToken();
					symbol = st.nextToken();
					locustag = st.nextToken();
					nick = st.nextToken();
					dbxrefs = st.nextToken();
					dbxrefs = dbxrefs.trim();
					chr = st.nextToken();
					map = st.nextToken();
					desc = st.nextToken();
					type = st.nextToken();
					xdb = "-";
					xdbid = "-";
					tax_final = "";
					tax_final = mapTerm(taxid,tax1,tax2,lt);
					xdb = mapTerm(taxid,tax1,db,lt);
					xdb_final = "-";
					if (!geneid.equals("-") && !tax_final.equals(""))
					{
							if (!dbxrefs.equals("-"))
							{
								tmp1 = "";
								tmp2 = "";
								StringTokenizer st2 = new StringTokenizer(dbxrefs,"|:");
								while (st2.countTokens()>1)
								{
									tmp1 = st2.nextToken();
									tmp2 = st2.nextToken();
									if (tmp1.equalsIgnoreCase(xdb))
									{
										xdbid = tmp2;
										xdb_final = xdb;
									}
								}
							}
							if(!tax_final.equals(""))
							{
								count++;
								out.write(geneid+"\t"+tax_final+"\t"+symbol+"\t"+desc+"\t"+locustag+"\t"+xdbid+"\t"+xdb_final+"\t"+chr+"\t"+map+"\t"+type+"\n");
							}
					}
				}
			}
			System.out.println("count all genes:"+count_all);
			System.out.println("count selected genes:"+count);
			out.flush();
			out.close();
			b1.close();
		}
		catch(Exception e)
	  	{
	    	System.out.println("Exception of processing NCBI gene info.");
	    	System.out.println("error message:"+e.getMessage());
		}
	}
//// process geneid mapping files (3 methods)

    public static void processIDMappingMore(String InputFileList, String species, String gene_info, String uniprot2gene, String output)
    {//this method is to increase the uniprot2genneid mapping by mapping uniprot to species specific ID eg ZFIN gene id then use species specific ID to map to NCBI geneid based on gene information table, some species largely benefit from this extra mapping eg zebrafish, mouse and rat
		try
		{
			FileWriter out = new FileWriter(output);
            File f2 = new File(species);
            BufferedReader b2 = new BufferedReader(new FileReader(f2));
            String[] taxids1 = new String[100];
            String[] taxids2 = new String[100];
            String[] idtypes = new String[100];
            int count_species = fillStringArray(f2,taxids1,taxids2,idtypes);
            System.out.println("count species info table:"+count_species);
			String strr = "";
            File f1 = new File(InputFileList);
            BufferedReader b1 = new BufferedReader(new FileReader(f1));
            HashMap<String,String> protein2otherid = new HashMap<String,String>();
            int count=0;
            int count_other = 0;
            boolean need = false;
            String strr0 = "";
            String uniprot = "";
            String uniprot2 = "";
            String taxid = "";
            String gene = "";
            String type = "";
            String geneid = "";
            String geneid2 = "";
            String otherid = "";
            String info = "";
            int count_more =0;
            int count_id = 0;
            int count_unmapped = 0;
            int count_mapped = 0;
            String file_name = "";
            while( (strr=b1.readLine())!=null)//process uniprot ID mapping file and select uniprot to species specific IDs eg FlyBase, RGD...
            {
                if(strr.length()==0) continue;
                StringTokenizer st = new StringTokenizer(strr,"\t");
                if (st.countTokens()>1)
                {
					file_name = st.nextToken();
				}
                File f_result = new File(file_name);
                BufferedReader b4 = new BufferedReader(new FileReader(f_result));
                String str2 = "";
                while( (str2=b4.readLine())!=null)//read each file line by line
                {
					if (str2.length()==0 ) continue;
					count++;
					StringTokenizer st2 = new StringTokenizer(str2,"\t");
					uniprot= "";
					type = "";
					info = "";
					uniprot= st2.nextToken();
					type = st2.nextToken();
					info = st2.nextToken();
					need = false;
					need = searchTerm(type,idtypes,count_species);
					if (need)
					{
						if (isNumeric(info)) { info=type+":"+info;}
						protein2otherid.put(uniprot,info);
						count_other++;
					}
				}
				b4.close();
			}
            HashMap<String,String> gene2otherid = new HashMap<String,String>();//process gene info file to get geneid to species specific id
            File f3 = new File(gene_info);
            BufferedReader b3 = new BufferedReader(new FileReader(f3));
            String strr3 = "";
            int count2 = 0;
            while( (strr3=b3.readLine())!=null)//read the gene_info table then populate gene2otherid
            {
				otherid = "";
				type = "";
				String[] tmp = strr3.split("\t");
				gene = tmp[0];
				otherid = tmp[5];
				type = tmp[6];
				if(isNumeric(otherid))
				{
					otherid = type+":"+otherid;
					gene2otherid.put(otherid,gene);
					count2++;
				}
				else if (!otherid.equals("-"))
				{
					gene2otherid.put(otherid,gene);
					count2++;
				}
			}
			File f0 = new File(uniprot2gene);
			BufferedReader b0=new BufferedReader(new FileReader(f0));
            while( (strr0=b0.readLine())!=null)//
            {
				StringTokenizer st0 = new StringTokenizer(strr0,"\t");
				if (st0.countTokens()>2)
				{
					uniprot = st0.nextToken();
					taxid = st0.nextToken();
					geneid = st0.nextToken();
					geneid2="";
					count_id++;
				}
				else {System.out.println("error:uniprot2tax/geneid file");}
				if (geneid.equals("null") || geneid.equals(""))
				{
					count_unmapped++;
					StringTokenizer st1 = new StringTokenizer(uniprot,"-");
					uniprot2=st1.nextToken();
					otherid = protein2otherid.get(uniprot);
					if (otherid == null) {otherid = protein2otherid.get(uniprot2);}
					if (otherid!=null)
					{
						geneid2=gene2otherid.get(otherid);
					}
				}
				if (geneid2!=null && !geneid2.equals(""))
				{
					count_more++;
					out.write(uniprot+"\t"+taxid+"\t"+geneid2+"\n");
				}
				else { out.write(uniprot+"\t"+taxid+"\t"+geneid+"\n");}
			}
            out.flush();
            out.close();
             b0.close();b1.close();b3.close();
            System.out.println("count all uniprot ids: "+count_id);
            System.out.println("count unmapped: "+count_unmapped);
            System.out.println("count more uniprot2geneid mapped: "+count_more);
        }
        catch(Exception e)
        {
			System.out.println("Exception in the extra step of mapping uniprot to gene");
			System.out.println(e.getMessage());
		}
    }
    public static void processIDMapping(String InputFileList, String output)
    {//this method maps uniprot to NCBI geneid based on the geneid mapping files provided by UniProt
		try
		{
			Random r= new Random();
			String ts = r.nextInt(1000)+"_uniprot2species.txt";
			System.out.println("file_name: "+ts);
			UniProt2TaxID(InputFileList, ts);
            File f = new File(InputFileList);
            int count =0;
            BufferedReader b1 = new BufferedReader(new FileReader(f));
            BufferedReader b2 = null;
            String strr="";
            String strr0="";
            String file_name = "";
            String taxid = "";
            String uniprot= "";
            String uniprot2= "";
            String geneid= "";
            String type = "";
            String info = "";
            boolean hasGeneID = false;
            FileWriter out = new FileWriter(output);//big result file
            HashMap<String,String> protein2geneid = new HashMap<String,String>();
            while( (strr=b1.readLine())!=null)//read the file name one by one
            {
                if(strr.length()==0) continue;
                StringTokenizer st = new StringTokenizer(strr,"\t");
                if (st.countTokens()>1)
                {
					file_name = st.nextToken();
					taxid = st.nextToken();
				}
                count++;
                File f_result = new File(file_name);
                b2 = new BufferedReader(new FileReader(f_result));
                String str2 = "";
                while( (str2=b2.readLine())!=null)//read each file line by line
                {
					if (str2.length()==0 ) continue;
					count++;
					StringTokenizer st2 = new StringTokenizer(str2,"\t");
					uniprot= "";
					type = "";
					info = "";
					uniprot= st2.nextToken();
					type = st2.nextToken();
					info = st2.nextToken();
					if (type.equals("GeneID"))
					{
						if(!protein2geneid.containsKey(uniprot))
						{
							protein2geneid.put(uniprot,info);
						}
						else
						{
							geneid=protein2geneid.get(uniprot);
							geneid = geneid.concat(","+info);
							protein2geneid.put(uniprot,geneid);
						}
					}
				}
			}
			File f0 = new File(ts);
			BufferedReader b0=null;
            b0= new BufferedReader(new FileReader(f0));
            while( (strr0=b0.readLine())!=null)//read the file name one by one
            {
				StringTokenizer st0 = new StringTokenizer(strr0,"\t");
				uniprot2 = st0.nextToken();
				taxid = st0.nextToken();
				geneid = "";
				geneid = protein2geneid.get(uniprot2);
				if (geneid == null)//eg. "O15504-1" does not return geneid, check one more time with "O15504"
				{
					StringTokenizer st1 = new StringTokenizer(uniprot2,"-");
					geneid = protein2geneid.get(st1.nextToken());
				}

				out.write(uniprot2+"\t"+taxid+"\t"+geneid+"\n");
			}
            out.flush();
            out.close();
            b0.close();b1.close();b2.close();
            System.out.println("count: "+count);
        }
        catch(Exception e)
        {
			System.out.println("Exception extracting uniprot2geneid mapping using uniprot id mapping files");
			System.out.println(e.getMessage());
		}
    }

    public static void UniProt2TaxID(String InputFileList, String output)
    {
		try
		{
            File f = new File(InputFileList);
            int count =0;
            BufferedReader b1=new BufferedReader(new FileReader(f));
            BufferedReader b2=null;
            String strr="";
            String file_name = "";
            String taxid = "";
            String uniprot= "";
            FileWriter out = new FileWriter(output);
            HashMap<String,String> protein2taxid = new HashMap<String,String>();
            while( (strr=b1.readLine())!=null)//read the file name one by one
            {
                if(strr.length()==0) continue;
                StringTokenizer st = new StringTokenizer(strr,"\t");
                if (st.countTokens()>1)
                {
					file_name = st.nextToken();
					taxid = st.nextToken();
				}
				System.out.println("file_name: "+file_name);
                count++;
                File f_result = new File(file_name);
                b2 = new BufferedReader(new FileReader(f_result));
                String str2 = "";
                while( (str2=b2.readLine())!=null)//read eah file line by line
                {
					if (str2.length()==0) continue;
					count++;
					StringTokenizer st2 = new StringTokenizer(str2,"\t");
					uniprot= "";
					uniprot= st2.nextToken();
					protein2taxid.put(uniprot,taxid);
				}
			}
			String[] keys = protein2taxid.keySet().toArray(new String[0]);
			System.out.println(keys.length);
			for  (String key: keys){ out.write(key+"\t"+protein2taxid.get(key)+"\n");}
            out.flush();
            out.close();
            b1.close();b2.close();
            System.out.println("count: "+count);
        }
        catch(Exception e)
        {
			System.out.println("Exception in UniProt2TaxID");
			System.out.println(e.getMessage());
		}
    }

/////step 3, process QFO ortholog prediction files (3 methods), some files only contain protein pairs of uniprot ids while compara uses internal geneid
    public static void processQFO(String input, String ids, String output, int option)
    {//
 		try
		{
            File f = new File(input);
            BufferedReader b1 = new BufferedReader(new FileReader(f));
            String strr="";
            String filename = "";
            String method = "";
            int count = 0;
 			Random r= new Random();
			String tmp_file = r.nextInt(1000)+"_"+output;
            while( (strr=b1.readLine())!=null)//read the file line by line
            {
				if(strr.length()==0) continue;
                StringTokenizer st = new StringTokenizer(strr,"\t");
                if(st.countTokens()>1)
                {
					count++;
					method = st.nextToken();
					filename = st.nextToken();
					if (method.equalsIgnoreCase("Compara"))
					{
						processCompara(filename,tmp_file);
					}
					else
					{
						processPair(filename,method,tmp_file);
					}
				}
			}
            filename = filterFile(tmp_file,ids, option);
            removeRedundancy(filename,output);
            System.out.println("option: "+option);
            System.out.println("count all: "+count);
            b1.close();
        }
        catch(Exception e)
        {
			System.out.println("Exceptions of parsing QFO file");
			System.out.println(e.getMessage());
		}
    }
    public static void processPair(String Input, String method, String Output)
    {//majority QFO files provides uniprot pairs in the first two columns and this function is just to order the IDs (eg. synchronize A-B,B-A to A-B) for integration
 		try
		{
            File f = new File(Input);
            BufferedReader b1= new BufferedReader(new FileReader(f));
            String strr="";
            String proteinid1 = "";
            String proteinid2 = "";
            int count = 0;
            FileWriter out = new FileWriter(Output,true);
            while( (strr=b1.readLine())!=null)//read the file line by line
            {
				if(strr.length()==0) continue;
                StringTokenizer st = new StringTokenizer(strr,"\t");
                if(st.countTokens()>1)
                {
					count++;
					proteinid1 = st.nextToken();
					proteinid2 = st.nextToken();
					if (proteinid1 == null || proteinid2 == null) {continue;}
                	if (proteinid1.compareTo(proteinid2)>0)
                	{
						out.write(proteinid1+"\t"+proteinid2+"\t"+method+"\n");
					}
					else if (proteinid1.compareTo(proteinid2)<0)
					{
						out.write(proteinid2+"\t"+proteinid1+"\t"+method+"\n");
					}
					else {continue;}
				}
			}
            out.flush();
            out.close();
            b1.close();
            System.out.println("method: "+method);
            System.out.println("count all: "+count);
        }
        catch(Exception e)
        {
			System.out.println("Exceptions parsing pair file");
			System.out.println(e.getMessage());
		}
    }

    public static void processCompara(String Input, String Output)
    {//process Compara xml file
		try
		{
            File f = new File(Input);
            int count =0;
            BufferedReader b1 = new BufferedReader(new FileReader(f));
            String strr="";
            String name = "";
            FileWriter out = new FileWriter(Output,true);
            HashMap<String,String> gene2protein = new HashMap<String,String>();
            String[] tmp = new String[10];
            String[] tmp2 = new String[20];
            int count_protein = 0;
            int count_pair = 0;
            String geneid = "";
            String geneid1 = "";
            String geneid2 = "";
            String proteinid = "";
            String proteinid1 = "";
            String proteinid2 = "";
            while( (strr=b1.readLine())!=null)//read the file name one by one
            {
                if(strr.length()==0) continue;
                count++;
                if (strr.indexOf("gene id")>0 || strr.indexOf("protId")>0)
                {
					tmp=strr.split("\"");
                	geneid = tmp[1];
                	proteinid= tmp[3];
                	count_protein++;
                	gene2protein.put(geneid,proteinid);
				}
				else if (strr.indexOf("orthologGroup id")>0)
				{
					count_pair++;
					tmp2=strr.split("\"");
					geneid1="";
					geneid1="";
	               	geneid1 = tmp2[7];
                	geneid2 = tmp2[9];
                	proteinid1= "";
                	proteinid2= "";
                	proteinid1=gene2protein.get(geneid1);
                	proteinid2=gene2protein.get(geneid2);
                	if (proteinid1.compareTo(proteinid2)>0)
                	{
						out.write(proteinid1+"\t"+proteinid2+"\tCompara\n");
					}
					else if (proteinid1.compareTo(proteinid2)<0)
					{
						out.write(proteinid2+"\t"+proteinid1+"\tCompara\n");
					}
				}
			}
            out.flush();
            out.close();
            b1.close();
            System.out.println("method: Compara");
            System.out.println("count_protein: "+count_protein);
            System.out.println("count_pair: "+count_pair);
        }
        catch(Exception e)
        {
			System.out.println("Exceptions of Compara parser");
			System.out.println(e.getMessage());
		}
    }// end of pasting method

/////step 4 (optinal): process ortholog predictions from HGNC and ZFIN
	static void processMore(String inputfile_list, String gene_info, String output)
	{
		try
		{
			HashMap<String,String> gene2tax = new HashMap<String,String>();
			HashMap<String,String> otherid2gene = new HashMap<String,String>();
			File f = new File(inputfile_list);
			BufferedReader b = new BufferedReader(new FileReader(f));
			BufferedReader b0 = null;
			BufferedReader b2 = null;
			String strr0 = "";
			FileWriter out = new FileWriter(output,true);
			String strr = "";
			int count_selected = 0;
			int count = 0;
			String method = "";
			String filename = "";
			int col = 0;
			String idtype = "";
			String id1 = "";
			String id2 = "";
			String tax1 = "";
			String tax2 = "";
			String geneid1 = "";
			String geneid2 = "";
			String strr2 = "";
			String taxs = "";
			String tax = "";
			String geneid = "";
			String otherid = "";
			String otheridtype = "";
			String key = "";
	        while (((strr=b.readLine())!=null))
	        {
				System.out.println(strr);
				count=0;
				count_selected=0;
				StringTokenizer st = new StringTokenizer(strr,"\t");
				if (st.countTokens()>4)
				{
					method = st.nextToken();
					filename = st.nextToken();
					col = Integer.parseInt(st.nextToken());
					idtype = st.nextToken();
					taxs = st.nextToken();
					StringTokenizer st0 = new StringTokenizer(taxs,"-");
					if (st0.countTokens()>1)
					{
						tax1 = st0.nextToken();
						tax2 = st0.nextToken();
					}
				}
				else
				{
					System.out.println("failed, missing info in the annotation file");
					return;
				}
				File f0 = new File(gene_info);
				b0 = new BufferedReader(new FileReader(f0));
	        	while (((strr0=b0.readLine())!=null))//populate the hashmaps based on gene_info file
	        	{
					count++;
					String[] tmp_gene = strr0.split("\t");
					tax = tmp_gene[1];
					geneid = tmp_gene[0];
					otherid= tmp_gene[5];
					otheridtype = tmp_gene[6];
					if (tax.equals(tax1)||tax.equals(tax2) )
					{

						gene2tax.put(geneid,tax);
						if (isNumeric(otherid)&& !otherid.equals("-"))
						{
							otherid = otheridtype+":"+otherid;
							otherid2gene.put(otherid,geneid);
							count_selected++;
						}
						else if (!otherid.equals("-"))
						{
							otherid2gene.put(otherid,geneid);
							count_selected++;
						}
					}
				}
				System.out.println("count gene info file:"+count);
				System.out.println("count selected:"+count_selected);
				count=0;
				count_selected=0;
				File f2 = new File(filename);
				b2 = new BufferedReader(new FileReader(f2));
				while (((strr2=b2.readLine())!=null))
				{
					count++;
					String[] tmp = strr2.split("\t");
					id1 = tmp[0];
					if (tmp.length<col) {continue;}
					id2 = tmp[col-1];
					if ((method.equals("HGNC") || method.equals("RGD"))&& !id2.equals(""))
					{
						geneid1 = otherid2gene.get(id1);
						StringTokenizer st2 = new StringTokenizer(id2,",");
						while (st2.hasMoreTokens())
						{
							geneid2 = otherid2gene.get(st2.nextToken());
							tax1= gene2tax.get(geneid1);
							tax2= gene2tax.get(geneid2);
							if (geneid1!=null && geneid2 != null && tax1 != null && tax2 != null)
							{
								out.write(tax1+"\t"+geneid1+"\t"+tax2+"\t"+geneid2+"\t"+method+"\n");
								out.write(tax2+"\t"+geneid2+"\t"+tax1+"\t"+geneid1+"\t"+method+"\n");
								count_selected++;
							}
						}
					}
					else if (method.equals("ZFIN")&& !id2.equals(""))
					{
						geneid1 = otherid2gene.get(id1);
						tax1 = gene2tax.get(geneid1);
						if (idtype.equals("geneid"))
						{
							geneid2= id2;
						}
						else
						{
							geneid2 = otherid2gene.get(id2);
						}
						tax2=gene2tax.get(geneid2);
						if (geneid1!=null && geneid2 != null && tax1 != null && tax2 != null)
						{
							out.write(tax1+"\t"+geneid1+"\t"+tax2+"\t"+geneid2+"\t"+method+"\n");
							out.write(tax2+"\t"+geneid2+"\t"+tax1+"\t"+geneid1+"\t"+method+"\n");
							count_selected++;
						}
					}
				}
				System.out.println("all records:"+count);
				System.out.println("selected records:"+count_selected);
			}
			out.flush();
			out.close();
			b0.close();b2.close();b.close();
		}
		catch(Exception e)
		{
			System.out.println("Exception processing more pairs from HGNC/ZFIN");
			System.out.println(e.getMessage());
		}
	}



/////small utility 1 - combine files
    public static void combine(String InputFileList, String OutputFile)
    {
		try
		{

            File f = new File(InputFileList);
            int count =0;
            BufferedReader b1=null;
            b1= new BufferedReader(new FileReader(f));
            String strr="";
            FileWriter out = new FileWriter(OutputFile);//combined result file
            while( (strr=b1.readLine())!=null)//read the file name one by one
            {
                if(strr.length()==0) continue;
                System.out.println("file_name: "+strr);
                File f_result = new File(strr);
                count++;
                BufferedReader b2=null;
                b2 = new BufferedReader(new FileReader(f_result));
                String str2 = "";
                while( (str2=b2.readLine())!=null)
                {
					if (str2.length()==0) continue;
					out.write(str2+"\n");
				}
            }// end of the while loop reading the names of the medline files
            out.flush();
            out.close();
            b1.close();
            System.out.println("count: "+count);
        }catch(Exception e){System.out.println("Exception combining files!!!");}
    }
/////small utility 2 - remove redundancy
	static void removeRedundancy(String inputfile, String output)
	{
		try
		{
			HashSet hs = new HashSet();
			File f = new File(inputfile);
			BufferedReader b = new BufferedReader(new FileReader(f));
			FileWriter out = new FileWriter(output);
			String strr = "";
			int count = 0;
			int count_selected = 0;
			boolean want = false;
	        while (((strr=b.readLine())!=null))
	        {
				count++;
				want = hs.contains(strr);
				if (!want)
				{
					hs.add(strr);
					out.write(strr+"\n");
					count_selected++;
				}
			}
			b.close();
			out.flush();
			out.close();
			System.out.println("all records:"+count);
			System.out.println("selected records:"+count_selected);
		}
		catch(Exception e)
		{
			System.out.println("Exception removing redundancy");
			System.out.println(e.getMessage());
		}
	}

/////small utility 3 - filter files
    public static String filterFile(String input, String ids, int option)
    {//filter the uniprot pairs based on the selected species and with the option to map to geneid (option=1).  if option<>1, uniprotids will not be mapped to geneid
		Random r= new Random();
		String tmp_file2 = r.nextInt(1000)+"_"+input;
		try
		{
            File f1 = new File(input);
            int count =0;
            int count1 =0;
            BufferedReader b1 = new BufferedReader(new FileReader(f1));
            File f2 = new File(ids);
            BufferedReader b2 = new BufferedReader(new FileReader(f2));
            String strr1="";
            String strr2="";
            String uniprot = "";
            String uniprot1 = "";
            String uniprot2 = "";
            String method = "";
            String geneid = "";
            String geneid1 = "";
            String geneid2 = "";
            String gene1 = "";
            String gene2 = "";
            String tax = "";
            String tax1 = "";
            String tax2 = "";
            FileWriter out = new FileWriter(tmp_file2);
            HashMap<String,String> protein2taxid = new HashMap<String,String>();
            HashMap<String,String> protein2geneid = new HashMap<String,String>();
            while( (strr2=b2.readLine())!=null)//read the file of ids
            {
				if (strr2.length()==0) continue;
          		StringTokenizer st2 = new StringTokenizer(strr2,"\t");
                if (st2.countTokens()>2)
                {
					uniprot = st2.nextToken();
					tax = st2.nextToken();
					geneid = st2.nextToken();
					protein2taxid.put(uniprot,tax);
					count++;
					if (!geneid.equals("null")) {protein2geneid.put(uniprot,geneid);count1++;}
				}
			}
			System.out.println("ID file:"+count);
			System.out.println("ID file with geneid:"+count1);
            while( (strr1=b1.readLine())!=null)//read the file of uniprot pairs
            {
				if (strr1.length()==0) continue;
          		StringTokenizer st1 = new StringTokenizer(strr1,"\t");
                if (st1.countTokens()>2)
                {
					uniprot1 = st1.nextToken();
					uniprot2 = st1.nextToken();
					method = st1.nextToken();
					tax1= protein2taxid.get(uniprot1);
					tax2= protein2taxid.get(uniprot2);
					if (option==1 && tax1!=null && tax2!=null)
					{
						geneid1= protein2geneid.get(uniprot1);
						geneid2= protein2geneid.get(uniprot2);
						if (geneid1!=null && geneid2!=null)
						{
							StringTokenizer stt1 = new StringTokenizer(geneid1,",");
							while (stt1.hasMoreTokens())
							{
								gene1=stt1.nextToken();
								StringTokenizer stt2 = new StringTokenizer(geneid2,",");
								while (stt2.hasMoreTokens())
								{
									gene2 = stt2.nextToken();
									out.write(tax1+"\t"+gene1+"\t"+tax2+"\t"+gene2+"\t"+method+"\n");
									out.write(tax2+"\t"+gene2+"\t"+tax1+"\t"+gene1+"\t"+method+"\n");
								}
							}
						}
					}
					else if (tax1!=null && tax2!=null)
					{
						out.write(tax1+"\t"+uniprot1+"\t"+tax2+"\t"+uniprot2+"\t"+method+"\n");
						out.write(tax2+"\t"+uniprot2+"\t"+tax1+"\t"+uniprot1+"\t"+method+"\n");
					}
				}
			}
			b1.close();
			b2.close();
			out.close();
		}
		catch(Exception e)
	  	{
	    	System.out.println("exception of filtering UniProt pairs based on selected species");
	    	System.out.println("error message:"+e.getMessage());
		}
		return (tmp_file2);
	}

/////other small utility
    public static boolean isNumeric(String str)
    {
	   try
	   {
	     Double.parseDouble(str);
	     return true;
	   }
	   catch(NumberFormatException e)
	   {
	     return false;
	   }
	}
/////other small utility
	static String mapTerm(String term, String[] term1, String[] term2, int len)
	{
		String tmp = "";
		String ret = "";
		boolean result = false;
		for (int t = 0; t<len; t++)
		{
			tmp = term1[t];
			if (tmp.equals(term))
			{
				result = true;
				ret = term2[t];
				break;
			}
		}
		return ret;
	}
/////other small utility
	static boolean searchTerm(String term, String[] terms, int len)
	{
		String tmp = "";
		boolean ret = false;
		for (int t = 0; t<len; t++)
		{
			tmp = terms[t];
			if (tmp.equals(term))
			{
				ret = true;
				break;
			}
		}
		return ret;
	}
/////other small utility
    public static int fillStringArray(File f,String[] arr1, String[] arr2, String[] arr3)
    {
        int count = 0;
        String strr="";
        BufferedReader b1=null;
        try
        {
            b1= new BufferedReader(new FileReader(f));
            while( (strr=b1.readLine())!=null)
            {
				StringTokenizer st = new StringTokenizer(strr,"\t");
				if (st.countTokens()==1)
				{
                	arr1[count]=st.nextToken();
                	count++;
				}
 				if (st.countTokens()==2)
 				{
                 	arr1[count]=st.nextToken();
                 	arr2[count]=st.nextToken();
                 	count++;
				}
 				if (st.countTokens()==3)
 				{
                 	arr1[count]=st.nextToken();
                 	arr2[count]=st.nextToken();
                 	arr3[count]=st.nextToken();
                 	count++;
				}
			}
        }
        catch(Exception e)
        {
            System.out.println("Exceptions in filling string arrays.");
        }
        return count;
    }


/////step 5: count vote and make ortholog pair best table and annotate the table (2 methods)
     public static void makeOrthologPairBest(String InputFile, String OutputFile)
     {
 		try
 		{
			HashMap<String,Integer> pair = new HashMap<String, Integer>(); //this is hashmap of unique gene pairs
 			String t1 = "";
 			String t2 = "";
 			String g1 = "";
 			String g2 = "";
 			String key_string = "";
 			int value = 0;
            File f = new File(InputFile);
            int count =0;
            BufferedReader b1 = new BufferedReader(new FileReader(f));
            String strr="";
            boolean result;
            FileWriter out = new FileWriter(OutputFile);
            while( (strr=b1.readLine())!=null)//read the file name one by one
            {
				if(strr.length()==0) continue;
                count++;
                key_string = "";
                StringTokenizer st = new StringTokenizer(strr,"\t");
 				if(st.countTokens()>3)
 				{
 					t1 = st.nextToken();
 					g1 = st.nextToken();
 					t2 = st.nextToken();
 					g2 = st.nextToken();
 					key_string = t1+"-"+g1+"-"+t2+"-"+g2;
 					result = false;
 					if (count==1)
 					{
 						pair.put(key_string,1);
 					}
 					else
 					{
 						result = pair.containsKey(key_string);
 						if (!result)
 						{
 							pair.put(key_string,1);
 						}
 						else
 						{
 							value = pair.get(key_string);
 							pair.remove(key_string);
 							value++;
 							pair.put(key_string,value);
 						}
 					}
 				}
 			}
 			for (String key : pair.keySet())
 			{

 				StringTokenizer st2 = new StringTokenizer(key,"-");
 				while (st2.hasMoreTokens())
 				{
 					out.write(st2.nextToken()+"\t");
 				}
 				out.write(pair.get(key)+"\n");
 			}
 			b1.close();
            out.flush();
            out.close();
         }catch(Exception e){System.out.println("Exception of counting the vote");}
    }
    public static void annotateOrthologPairBest(String InputFile, String OutputFile)
    {
		String strr="";
		try
		{
			HashMap<String,Integer> pair = new HashMap<String, Integer>(); //this is hashmap of unique gene pairs
			String t1 = "";
			String t2 = "";
			String g1 = "";
			String g2 = "";
			String key_string = "";
			String key_string1 = "";
			String key_string2 = "";
			int value = 0;
			int score = 0;
			int value1 = 0;
			int value2 = 0;
			String best1 = "";
			String best2 = "";
			String rank = "";
            File f = new File(InputFile);
            int count =0;
            int count0 =0;
            int count_pair =0;
            BufferedReader b1 = new BufferedReader(new FileReader(f));
            BufferedReader b2 = new BufferedReader(new FileReader(f));
            boolean result;
            FileWriter out = new FileWriter(OutputFile);
            FileWriter out2 = new FileWriter("hashmap_pair.txt");
            while( (strr=b1.readLine())!=null)
            {
                if(strr.length()==0) continue;
                count++;
                key_string = "";
				StringTokenizer st = new StringTokenizer(strr,"\t");
				if(st.countTokens()>4)
				{
					t1 = st.nextToken();
					g1 = st.nextToken();
					t2 = st.nextToken();
					g2 = st.nextToken();
					score = Integer.parseInt(st.nextToken());
					key_string = t1+"-"+g1+"-"+t2;
					if (count==1)
					{
						pair.put(key_string,score);
						count_pair++;
					}
					else
					{
						result = pair.containsKey(key_string);
						if (!result)
						{
							pair.put(key_string,score);
							count_pair++;
						}
						else
						{
							value = pair.get(key_string);
							if (score > value)
							{
								pair.put(key_string,score);
							}
						}
					}
				}
			}
			System.out.println("populate best score hashmap:"+count_pair);
			out2.write("hashmap:"+pair);
			out2.close();
            while( (strr=b2.readLine())!=null)//read the file name one by one
            {
                if(strr.length()==0) continue;
                count0++;
                key_string = "";
				StringTokenizer st = new StringTokenizer(strr,"\t");
				if(st.countTokens()>4)
				{
					t1 = st.nextToken();
					g1 = st.nextToken();
					t2 = st.nextToken();
					g2 = st.nextToken();
					value1 = 0;
					value2 = 0;
					key_string1 = "";
					key_string2 = "";
					score = Integer.parseInt(st.nextToken());
					key_string1 = t1+"-"+g1+"-"+t2;
					value1 = pair.get(key_string1);
					key_string2 = t2+"-"+g2+"-"+t1;
					value2 = pair.get(key_string2);
					best1 = "No";
					best2 = "No";
					rank = "low";
					if (score==value1)
					{
						best1 = "Yes";
					}
					if (score == value2)
					{
						best2 = "Yes";
					}
					if (best1.equals("Yes") && best2.equals("Yes") && score>1)
					{
						rank = "high";
					}
					else if (((best1.equals("Yes") || best2.equals("Yes")) && score>1) || score>3)
					{
						rank = "moderate";
					}
					out.write(count0+"\t"+t1+"\t"+g1+"\t"+t2+"\t"+g2+"\t"+score+"\t"+best1+"\t"+best2+"\t"+rank+"\n");
				}
			}
			out.close();
			System.out.println("count:"+count);
			System.out.println("count1:"+count0);
        }
        catch(Exception e)
        {
			System.out.println("Exception when annotate the Ortholog_Pair table for rank and evaluation of best score");
			System.out.println(e.getMessage());
			System.out.println(strr);
		}
    }
    private static final int DEFAULT_ARRAY_SIZE = 12000000;
    public static void main(String[] args)
    {
        if (args.length < 1) {
            printUsage();
            System.exit(1);
        }

        String operation = args[0];

        try {
            switch (operation.toLowerCase()) {
                case "assemblegeneinfo":
                    handleAssembleGeneInfo(args);
                    break;
                case "processidmapping":
                    handleProcessIDMapping(args);
                    break;
                case "processidmappingmore":
                    handleProcessIDMappingMore(args);
                    break;
                case "processqfo":
                    handleProcessQFO(args);
                    break;
                case "processmore":
                    handleProcessMore(args);
                    break;
                case "removeredundancy":
                    handleRemoveRedundancy(args);
                    break;
                case "combine":
                    handleCombine(args);
                    break;
                case "makeorthologpairbest":
                    handleMakeOrthologPairBest(args);
                    break;
                case "annotateorthologpairbest":
                    handleAnnotateOrthologPairBest(args);
                    break;
                default:
                    System.err.println("Unknown operation: " + operation);
                    printUsage();
                    System.exit(1);
            }
        } catch (Exception e) {
            System.err.println("Error executing operation: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private static void printUsage() {
        System.out.println("Usage: java DIOPT_QFO <operation> [options]");
        System.out.println("\nOperations:");
        System.out.println("  assembleGeneInfo <ncbi gene inputfile> <species_file> <outputfile-gene_information>");
        System.out.println("  processIDMapping <list of ID mapping files from UniProt> <outputfile-uniprot2geneid>");
        System.out.println("  processIDMappingMore <list of ID mapping files from UniProt> <species_file> <gene_information> <current uniprot2gene> <output-uniprot2gene with more mapping>");
        System.out.println("  processQFO <list of selected method> <uniprot2gene with more mapping> <output-ortholog_pair> [option]");
        System.out.println("  processMore <list of more resources> <gene_information> <output-ortholog pair more>");
        System.out.println("  removeRedundancy <output-ortholog pair more> <output-ortholog pair more unique records>");
        System.out.println("  combine <list of ortholog pair files> <output-final ortholog pair>");
        System.out.println("  makeOrthologPairBest <final ortholog pair> <ortholog pair best>");
        System.out.println("  annotateOrthologPairBest <ortholog pair best> <final ortholog pair best with annotation>");
    }


    // Handler methods
    private static void handleAssembleGeneInfo(String[] args) throws IOException {
        if (args.length < 4) {
            System.err.println("Usage: assemblegeneinfo <ncbi gene inputfile> <species_file> <outputfile-gene_information>");
            System.exit(1);
        }
        assembleGeneInfo(args[1], args[2], args[3]);
    }

    private static void handleProcessIDMapping(String[] args) throws IOException {
        if (args.length < 3) {
            System.err.println("Usage: processIDMapping <list of ID mapping files from UniProt> <outputfile-uniprot2geneid>");
            System.exit(1);
        }
        processIDMapping(args[1], args[2]);
    }

    private static void handleProcessIDMappingMore(String[] args) throws IOException {
        if (args.length <6 ) {
            System.err.println("Usage: processIDMappingMore <list of ID mapping files from UniProt> <species_file> <gene_information> <current uniprot2gene> <output-uniprot2gene with more mapping>");
            System.exit(1);
        }
        processIDMappingMore(args[1], args[2], args[3], args[4],args[5]);
    }

    private static void handleProcessQFO(String[] args) throws IOException {
        if (args.length < 5) {
            System.err.println("Usage: processQFO <list of selected method> <uniprot2gene with more mapping> <output-ortholog_pair> [option]");
            System.exit(1);
        }
        processQFO(args[1], args[2],args[3],Integer.parseInt(args[4]));
    }

    private static void handleProcessMore(String[] args) throws IOException {
        if (args.length < 3) {
            System.err.println("processMore <list of more resources> <gene_information> <output-ortholog pair more>");
            System.exit(1);
        }
        processMore(args[1], args[2],args[3]);
	}

    private static void handleRemoveRedundancy(String[] args) throws IOException {
        if (args.length < 3) {
            System.err.println("Usage: removeRedundancy <output-ortholog pair more> <output-ortholog pair more unique records>");
            System.exit(1);
        }
        removeRedundancy(args[1], args[2]);
	}
    private static void handleCombine(String[] args) throws IOException {
        if (args.length < 3) {
            System.err.println("Usage: combine <list of ortholog pair files> <output-final ortholog pair>");
            System.exit(1);
        }
        combine(args[1], args[2]);
	}
    private static void handleMakeOrthologPairBest(String[] args) throws IOException {
        if (args.length < 3) {
            System.err.println("Usage: makeOrthologPairBest <final ortholog pair> <ortholog pair best>");
            System.exit(1);
        }
        makeOrthologPairBest(args[1], args[2]);
    }
    private static void handleAnnotateOrthologPairBest(String[] args) throws IOException {
        if (args.length < 3) {
            System.err.println("Usage: annotateOrthologPairBest <ortholog pair best> <final ortholog pair best with annotation>");
            System.exit(1);
        }
        annotateOrthologPairBest(args[1], args[2]);
    }

}