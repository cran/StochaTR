// Here we describe the main agents which realize the algorithm
// P.Bazovkin, 2010.

//#include <windows.h>		// Header File For Windows
#include <time.h>
//#include <GL/gl.h>			// Header File For The OpenGL32 Library
//#include <GL/glu.h>			// Header File For The GLu32 Library
//#include <GL/glext.h>		// Header File For The Glaux Library

// Routine with the static elements
#define INCL_OBJECTS
#include "StochaTRObjects.cpp"		// Header File For my private objects Library
//#undef INCL_OBJECTS

// General interface for algorithmic agents
class Agent
{
};


// ********************** Class for displaying the result *****************************************

	bool	keys[256];		// Array Used For The Keyboard Routine
	bool	active;			// Window Active Flag Set To TRUE By Default


class ResultAg : public Agent
{
public:

	list<Facet> trimmed_region;
	vector<Point> data_cloud;

	Point coord_center;         // Coordinates of the centroid

	void PrintSolution(vector<float> _optfac)
	{
		ofstream os("slp-WMTR.dat");

                // Special cases
                if(_optfac[0] == -11)
                {
                    os << "The solution is unlimited";
                    os.close();
                    return;
                }
                else if(_optfac[0] == -22)
                {
                    os << "No solution";
                    os.close();
                    return;
                }
		
		os << "The optimal solution of the SLP is: \n";
		
		os << "( ";
		for(int i = 0; i < _optfac.size(); i++)
		{
			os << _optfac[i] << ", ";
		}		
		os << " )\n";

		os.close();
	}


	ResultAg()
	{
	}


	
private:
	// technical variables
	int flborder;
	bool markfirst;

public:

	void GLSceneToExport()									// Here's Where We Do All The Drawing
	{
		// Only for d=3
		if(this->coord_center.coord.size() != 3) return;

		int currcol[3];

		//// ********************* Display points from the data cloud **************************

		
		ofstream osvisd( "tmp_vis_dt.dat");

		currcol[0] = 255;
		currcol[1] = 0;
		currcol[2] = 0;

		vector<Point>::iterator pointit;

		float fat = 0.1f;

		// Typing coordinates of the mean at the first place
		//osvisd << this->coord_center.coord[0] << " " << this->coord_center.coord[1] << " " << this->coord_center.coord[2] << endl;
		osvisd << 0 << " " << 0 << " " << 0 << endl;

		for(pointit = this->data_cloud.begin(); pointit != this->data_cloud.end(); pointit++)
		{
			// Displaying as spheres	
			osvisd << pointit->coord[0] << " " << pointit->coord[1] << " " << pointit->coord[2] << endl;

			// Displaying as pyramides	
			/* osvisd << pointit->coord[0] + fat << " " << pointit->coord[1] << " " << pointit->coord[2] << endl;					
				osvisd << pointit->coord[0] - fat << " " << pointit->coord[1] + fat << " " << pointit->coord[2] << endl;					
				osvisd << pointit->coord[0] - fat << " " << pointit->coord[1] - fat << " " << pointit->coord[2] - fat << endl;					

				osvisd << pointit->coord[0] + fat << " " << pointit->coord[1] << " " << pointit->coord[2] << endl;					
				osvisd << pointit->coord[0] - fat << " " << pointit->coord[1] + fat << " " << pointit->coord[2] << endl;					
				osvisd << pointit->coord[0] - fat << " " << pointit->coord[1] - fat << " " << pointit->coord[2] + fat << endl;					
				
				osvisd << pointit->coord[0] + fat << " " << pointit->coord[1] << " " << pointit->coord[2] << endl;					
				osvisd << pointit->coord[0] - fat << " " << pointit->coord[1] - fat << " " << pointit->coord[2] - fat << endl;					
				osvisd << pointit->coord[0] - fat << " " << pointit->coord[1] - fat << " " << pointit->coord[2] + fat << endl;					

				osvisd << pointit->coord[0] - fat << " " << pointit->coord[1] + fat << " " << pointit->coord[2] << endl;					
				osvisd << pointit->coord[0] - fat << " " << pointit->coord[1] - fat << " " << pointit->coord[2] - fat << endl;					
				osvisd << pointit->coord[0] - fat << " " << pointit->coord[1] - fat << " " << pointit->coord[2] + fat << endl; */					
		}

		osvisd.close();


		// ***********************************************************************************


		// ********************* Display the trimmed region **********************************

		ofstream osvisf("tmp_vis_fc.dat");

		list<Facet>::iterator facetit;

		for(facetit = this->trimmed_region.begin(); facetit != this->trimmed_region.end(); facetit++)
		{
			//if(facetit->normalvec.coord[0] < 0 || facetit->normalvec.coord[1] < 0 || facetit->normalvec.coord[2] < 0)
			//	continue;

			if( !facetit->truncated /*&&  ++stop<flborder*/  )
			{

		
				float peak = max(fabs(facetit->normalvec.coord[0]), max(fabs(facetit->normalvec.coord[1]), fabs(facetit->normalvec.coord[2])));
				float gr_sigma = 2;
				float gr_mu    = 0.5f;
				float colx = -facetit->normalvec.coord[0] / (gr_sigma * peak) + gr_mu;
				float coly = -facetit->normalvec.coord[1] / (gr_sigma * peak) + gr_mu;
				float colz = -facetit->normalvec.coord[2] / (gr_sigma * peak) + gr_mu;
				//glColor3f(colx, coly, colz);

				list<ExtremePoint>::iterator expointit;
				
				expointit = facetit->nodes.begin();
				osvisf << expointit->coord[0] << " " << expointit->coord[1] << " " << expointit->coord[2] << endl;		
				
				expointit++;
				osvisf << expointit->coord[0] << " " << expointit->coord[1] << " " << expointit->coord[2] << endl;					
				
				expointit++;
				osvisf << expointit->coord[0] << " " << expointit->coord[1] << " " << expointit->coord[2] << endl;
			}
		}

	

		for(facetit = this->trimmed_region.begin(); facetit != this->trimmed_region.end(); facetit++)
		{			
						
			//if(facetit->normalvec.coord[0] < 0 || facetit->normalvec.coord[1] < 0 || facetit->normalvec.coord[2] < 0)
			//	continue;


			if( facetit->truncated /*&& ++stop<flborder*/ )
			{

				float peak = max(fabs(facetit->normalvec.coord[0]), max(fabs(facetit->normalvec.coord[1]), fabs(facetit->normalvec.coord[2])));
				float gr_sigma = 2;
				float gr_mu    = 0.5f;
				float colx = -facetit->normalvec.coord[0] / (gr_sigma * peak) + gr_mu;
				float coly = -facetit->normalvec.coord[1] / (gr_sigma * peak) + gr_mu;
				float colz = -facetit->normalvec.coord[2] / (gr_sigma * peak) + gr_mu;

				list<ExtremePoint>::iterator expointit;
				list<ExtremePoint>::iterator expointitback;
				
				//glBegin(GL_POLYGON);

				//glColor3f(colx, coly, colz);

				vector<vector<float> > cvec(facetit->nodes.size());
				int j;
				for(j = 0, expointit = facetit->nodes.begin(); expointit != facetit->nodes.end(); expointit++, j++)
				{
					cvec[j] = expointit->coord;

				}

				// Hexagonal
				if(cvec.size() == 6)
				{
					osvisf << (cvec[0])[0] << " " << (cvec[0])[1] << " " << (cvec[0])[2] << endl;
					osvisf << (cvec[1])[0] << " " << (cvec[1])[1] << " " << (cvec[1])[2] << endl;
					osvisf << (cvec[2])[0] << " " << (cvec[2])[1] << " " << (cvec[2])[2] << endl;

					osvisf << (cvec[2])[0] << " " << (cvec[2])[1] << " " << (cvec[2])[2] << endl;
					osvisf << (cvec[3])[0] << " " << (cvec[3])[1] << " " << (cvec[3])[2] << endl;
					osvisf << (cvec[4])[0] << " " << (cvec[4])[1] << " " << (cvec[4])[2] << endl;

					osvisf << (cvec[0])[0] << " " << (cvec[0])[1] << " " << (cvec[0])[2] << endl;
					osvisf << (cvec[5])[0] << " " << (cvec[5])[1] << " " << (cvec[5])[2] << endl;
					osvisf << (cvec[4])[0] << " " << (cvec[4])[1] << " " << (cvec[4])[2] << endl;

					osvisf << (cvec[2])[0] << " " << (cvec[2])[1] << " " << (cvec[2])[2] << endl;
					osvisf << (cvec[0])[0] << " " << (cvec[0])[1] << " " << (cvec[0])[2] << endl;
					osvisf << (cvec[4])[0] << " " << (cvec[4])[1] << " " << (cvec[4])[2] << endl;
				}
				// Quad
				else
				{
					osvisf << (cvec[0])[0] << " " << (cvec[0])[1] << " " << (cvec[0])[2] << endl;
					osvisf << (cvec[1])[0] << " " << (cvec[1])[1] << " " << (cvec[1])[2] << endl;
					osvisf << (cvec[2])[0] << " " << (cvec[2])[1] << " " << (cvec[2])[2] << endl;

					osvisf << (cvec[2])[0] << " " << (cvec[2])[1] << " " << (cvec[2])[2] << endl;
					osvisf << (cvec[3])[0] << " " << (cvec[3])[1] << " " << (cvec[3])[2] << endl;
					osvisf << (cvec[0])[0] << " " << (cvec[0])[1] << " " << (cvec[0])[2] << endl;
				}

				
			}
		}

		osvisf.close();
		

		// ***********************************************************************************
		
		
		// ********************* Display edges of the trimmed region *************************

		ofstream osvisr("tmp_vis_rg.dat");
		
		//glColor3f(0.4f,0.4f,0.4f);						
		for(facetit = this->trimmed_region.begin(); facetit != this->trimmed_region.end(); facetit++)
		{			

			list<ExtremePoint>::iterator expointit, expointitn;
			expointit = facetit->nodes.begin();
			expointitn = expointit;
			expointitn++;
			for(; expointitn != facetit->nodes.end(); expointit++, expointitn++)
			{
				osvisr << expointit->coord[0] << " " << expointit->coord[1] << " " << expointit->coord[2] << endl;		
				osvisr << expointitn->coord[0] << " " << expointitn->coord[1] << " " << expointitn->coord[2] << endl;		
			}
			expointitn = facetit->nodes.begin();
			osvisr << expointit->coord[0] << " " << expointit->coord[1] << " " << expointit->coord[2] << endl;		
			osvisr << expointitn->coord[0] << " " << expointitn->coord[1] << " " << expointitn->coord[2] << endl;		
				
		}

		// A line from origin to centroid
		osvisr << 0 << " " << 0 << " " << 0 << endl;		
		osvisr << this->coord_center.coord[0] << " " << this->coord_center.coord[1] << " " << this->coord_center.coord[2] << endl;		


		osvisr.close();

		// ***********************************************************************************


	}



};

//*************************************************************************************************


// ********************** Class for building a trimmed region *************************************

class ProcessAg : public Agent
{
public:

	Permutation perm;             // Current permutation

	Point initial_center;         // Coordinates of the centroid

	string WMTD_type;              // Notion of a WMTR

protected:

	list<Facet> trimmed_region; // Trimmed region: set of the facets

	float alpha;                  // Zonoid depth

	list<Facet> queue;            // Queue of the not processed facets

	int   dim  ;                  // Dimension of the cloud
	int   num  ;                  // Number of points

	Vector p;                     // Current direction for support function

	int border_index;             // Index of the border point (affected by the depth)
	vector<float> weight;         // Weigts vector (non-decreasing function on the permutation)

	HashTable ohash_table;        // "Optimistic" hash table

	void NormalizeWeights()
	{
		// Smoothing too small weights
		for(int i = 0; i < num; i++)
		{
			if(this->weight[i] < 1.0f / (this->num * 100))
				this->weight[i] = 0.0f;
		}

		float sum_of_weights = 0;
		for(int i = 0; i < num; i++)
		{
			sum_of_weights += this->weight[i];
		}

		// If there is an incorrect alpha parametere
		if(sum_of_weights < 0.01f)
		{
			 this->weight[num-1] = 1;
			 sum_of_weights = 1;
		}

		for(int i = 0; i < num; i++)
		{
			this->weight[i] = this->weight[i] / sum_of_weights;
		}
	}

public:

	// Generates weight vector basing on the type of a WMTR
	void WeightsGenerator()
	{
		if(this->WMTD_type == "zonoid")
		{
			int zonbord = num - floor(alpha * num) - 1;
			this->border_index = zonbord;

			for(int i = 0; i < this->num; i++)
			{
				if( i<zonbord )
				{
					this->weight[i] = 0;
				}
				else if( i==zonbord)
				{
					this->weight[i] = (alpha * num - floor(alpha * num)) / (alpha * num);
				}
				else
				{
					this->weight[i] = 1.0f / (alpha * num);
				}
			}
		}
		else if(this->WMTD_type == "ECH") 
		{
			int beta = 1.0f / this->alpha;

			for(int j = num-1; j >= 0; j--)
			{
				if(j < beta-1)
				{
					this->weight[j]= 0.0f;
				}
				else if(j == num-1)
				{
					//this->weight[i]= 1.0f / Comb(num, beta);
					this->weight[j] = ((float)beta)/num;
				}
				else
				{
					this->weight[j]= weight[j+1] * (j+2-beta) / (j+1);
				}
			}
		}
		else if(this->WMTD_type == "ECH*") 
		{
			float beta = 1.0f / this->alpha;

			for(int i = 0; i < num; i++)
			{
				this->weight[i]= ( pow(i+1, beta) - pow(i, beta) ) / pow(num, beta);
			}
		}
		else if(this->WMTD_type == "geometrical") 
		{
			float gmult = (1-alpha) / (1-pow(alpha, num));

			for(int i = 0; i < num; i++)
			{
				this->weight[i]= pow(alpha, num-i-1) * gmult;
			}
		}
		else if(this->WMTD_type == "manual") 
		{
		}
		else
		{
			int zonbord = num - floor(alpha * num) - 1;
			this->border_index = zonbord;

			for(int i = 0; i < this->num; i++)
			{
				if( i<zonbord )
				{
					this->weight[i] = 0;
				}
				else if( i==zonbord)
				{
					this->weight[i] = (alpha * num - floor(alpha * num)) / (alpha * num);
				}
				else
				{
					this->weight[i] = 1.0f / (alpha * num);
				}
			}
		}

	}

	// Use only this constructor for this class
	ProcessAg(string _type, float _depth, int _dim, int _num, vector<Point> _cloud, vector<float> _manweights):
																p(_dim),
																initial_center(_dim),
																ohash_table(_num),
																weight(_num)
	{
		this->WMTD_type = _type;
		
		dim   = _dim;
		num   = _num;

		alpha = _depth;

		// ************** Normalization of the data ************

		vector<Point>::iterator it;

		// Calculate centroid
		for(it = _cloud.begin(); it != _cloud.end(); it++)
		{
			for(int j=0; j<dim; j++)
			{
				initial_center.coord[j] += it->coord[j];

			}

		}

		for(int j=0; j<dim; j++)
		{
			initial_center.coord[j] = initial_center.coord[j] / num;

		}

		//// Set centroid to 0
		//for(it = _cloud.begin(); it != _cloud.end(); it++)
		//{
		//	for(int j=0; j<dim; j++)
		//	{
		//		it->coord[j] = it->coord[j] - initial_center.coord[j];

		//	}

		//}
		// ************************************************************


		perm.points.assign(_cloud.begin(), _cloud.end());

		// Ordering points according to DEFAULT support vector
		perm.Support(p);

		if(this->WMTD_type == "general")
			weight = _manweights;
		else
			WeightsGenerator();

		NormalizeWeights();

		this->ohash_table.weight = this->weight;

	}

protected:

	// Working with the hash-table
	
	list<boost::dynamic_bitset<> > hash_table;
	
	static bool HashOrder(boost::dynamic_bitset<> bs1, boost::dynamic_bitset<> bs2)
	{
		//::ResumeCumulTime();
		// Warning: not optimal procedure
		int bsrange = bs1.size();
		for(int i = 0; i < bsrange; i++)
		{
			if(bs1[i] > bs2[i])
			{
				//::StopCumulTime();
				return true;
			}
			else if(bs1[i] < bs2[i])
			{
				//::StopCumulTime();
				return false;
			}
		}
		//::StopCumulTime();
		return true;

		/*int byterange = bs1.m_bits.size();
		for(int i = 0; i < byterange; i++)
		{
			if((int)bs1.m_bits[i] > (int)bs2.m_bits[i])
				return true;
			else if((int)bs1.m_bits[i] < (int)bs2.m_bits[i])
				return false;
		}
		return true;*/
	}
	
	static bool ReversedHashOrder(boost::dynamic_bitset<> bs1, boost::dynamic_bitset<> bs2)
	{
		return !HashOrder(bs1, bs2);
	}

	void MarkByHash(Facet _new_fac)
	{
		// Warning: not optimal insertion algorithm!
		list<boost::dynamic_bitset<> >::iterator hashit;

		bool inserted = false;
		
		boost::dynamic_bitset<> curr_hash = _new_fac.GetHashMap();

		if(this->hash_table.size() == 0)
		{
			this->hash_table.push_back(curr_hash);
		}
		else
		{
			for(hashit = this->hash_table.begin(); hashit != this->hash_table.end(); hashit++)
			{
				// Searching a place for insertion with saving the order
				if(ProcessAg::HashOrder(*hashit, curr_hash))
				{
					this->hash_table.insert(hashit, curr_hash);
					inserted = true;
					return;
				}
			}
			
			if(!inserted)
			{
				this->hash_table.insert(hash_table.end(), curr_hash);
			}
		}
	}

	// NOT key function (see ContainsFacet) 
	bool CheckInHash(Facet _new_fac)
	{
		return binary_search(this->hash_table.begin(), this->hash_table.end(), _new_fac.GetHashMap(), &ProcessAg::ReversedHashOrder);
	}

protected:
	
	// Factorial
	boost::uint64_t Fact(int _n)
	{
		if (_n <= 1)
		{
			return 1;
		}
		else
		{
			return _n*Fact(_n-1);
		}
	}

	// Calculates _n!/_b!
	boost::uint64_t FactShifted(int _n, int _b)
	{
		if (_n <= _b)
		{
			return 1;
		}
		else
		{
			return _n*FactShifted(_n-1, _b);
		}
	}

	// Number of combinations _n by _b
	// !Warning - large numbers!
	boost::uint64_t Comb(int _n, int _b)
	{
		return FactShifted(_n, max(_n-_b, _b)) / Fact(min(_n-_b, _b)); 
		//return Fact(_n)/(Fact(_n-_b) * Fact(_b)); 
	}

protected:

	// Argument: index permutation of extr point
	// Result  : coordinates of this extr point
	vector<float> CurrentExtremePoint(vector<int> _index_perm)
	{
		Point extreme(this->dim);

		for(int j = 0; j < this->num; j++)
		{
			for(int i = 0; i<this->dim; i++)
			{
				extreme.coord[i] += this->perm.points[_index_perm[j]].coord[i] * this->weight[j];
			}

		}

		return extreme.coord;


	}

	vector<float> CurrentExtremePoint_gener(vector<int> _index_perm, Permutation _perm)
	{
		Point extreme(this->dim);

		for(int j = 0; j < this->num; j++)
		{
			for(int i = 0; i<this->dim; i++)
			{
				extreme.coord[i] += _perm.points[_index_perm[j]].coord[i] * this->weight[j];
			}

		}

		return extreme.coord;

	}

private:

	// Returns whether current facet has already been included into the trimmed region 
	bool ContainsFacet(Facet _currfacet)
	{
		return binary_search(this->hash_table.begin(), this->hash_table.end(), _currfacet.GetHashMap(), &ProcessAg::ReversedHashOrder);
	}

protected:
public:

private:

	// Obtaining $d-1$ independent vectors from a set of $d$ points (in general position)
	vector<vector<float> > GenIndepVectors(list<int> _def_set)
	{
		vector<vector<float> > ind_bundle(this->dim); 

		list<int>::iterator iter;
		list<int>::iterator pivot_iter = _def_set.begin();

		int j = 0;
		for(iter = _def_set.begin(); iter != _def_set.end(); iter++)
		{
			if(iter != pivot_iter)
			{
				Vector vec(this->dim);
				vec.FromPointToPoint(perm.points[*pivot_iter], perm.points[*iter]);
				ind_bundle[j] = vec.coord;
				j++;
			}
		}

		return ind_bundle;
	}

	// Generate combination (curr_comb.count num) by increasing pos element index
	vector<int> NextCombination(vector<int> curr_comb, int n)
	{
		if(curr_comb.size() == 0) return curr_comb;

		vector<int> new_comb = curr_comb;
		int c = new_comb.size();
		
		int pos = c-1;
		bool found = false;
		do
		{
			if(new_comb[pos] < n-c+pos)
			{
				new_comb[pos]++;
				found = true;
				for(int j = pos+1; j < c; j++)
				{
					new_comb[j] = new_comb[pos] + j - pos;
				}
			}

			pos--;
		}while( pos >= 0 && found == false);

		if(!found)
		{
			new_comb[0] = -1;
		}

		return new_comb;
	}

protected:

	
	// Generate set of points
	list<int> NextSet(list<int> curr_comb, int n)
	{
		list<int>::iterator lit;
		lit = curr_comb.begin();

		vector<int> vcomb(dim);

		for(int u = 0; u < dim; u++)
		{
			vcomb[u] = *lit;
			lit++;
		}

		vcomb = NextCombination(vcomb, n);
		
		lit = curr_comb.begin();
		for(int u = 0; u < dim; u++)
		{
			*lit = vcomb[u];
			lit++;
		}

		return curr_comb;
	}


	void MarkByHashH(Facet _new_fac)
	{
		// Warning: not optimal insertion algorithm!
		list<boost::dynamic_bitset<> >::iterator hashit;

		bool inserted = false;
		
		boost::dynamic_bitset<> curr_hash = _new_fac.GetHashMapH(this->num);

		if(this->hash_table.size() == 0)
		{
			this->hash_table.push_back(curr_hash);
		}
		else
		{
			for(hashit = this->hash_table.begin(); hashit != this->hash_table.end(); hashit++)
			{
				// Searching a place for insertion with saving the order
				if(ProcessAg::HashOrder(*hashit, curr_hash))
				{
					this->hash_table.insert(hashit, curr_hash);
					inserted = true;
					//::StopCumulTime();
					return;
				}
			}
			
			if(!inserted)
			{
				this->hash_table.insert(hash_table.end(), curr_hash);
			}
		}
	}

	// Creates a matrix (d x d) with first K rows - defining vectors of _fac
	vector<vector<float> > GetDefVectors_kFace(kFace _fac)
	{
		vector<vector<float> > def_vecs(this->dim);
		int plz = 0;

		list<int>::iterator itset, itcard;
		for(itset = _fac.anchors.begin(), itcard = _fac.cardinals.begin();
			itset != _fac.anchors.end(); 
			itset++, itcard++)
		{
			for(int y = 0; y < *itcard - 1; y++)
			{
				Vector dvec(this->dim);
				dvec.FromPointToPoint(perm.points[_fac.index_perm[*itset + y]], perm.points[_fac.index_perm[*itset + y + 1]]);
				
				def_vecs[plz] = dvec.coord;
				plz++;
			}
		}
	
		return def_vecs;
	}

	// Returns whether current facet has already been included into the trimmed region 
	bool ContainsFacetH(Facet _currfacet)
	{
		bool fres = binary_search(this->hash_table.begin(), this->hash_table.end(), _currfacet.GetHashMapH(this->num), &ProcessAg::ReversedHashOrder);
		return fres;
	}

	list<DataVector> GetCriticalVectors(Facet _currfac)
	{
		// !!! Attention: _currfac.anchors should modified (with removed vector)

		// The first step - forming associating clusters

		list<int>     assoc_clust;            // Association clusters starting points in the pi()
		list<bool>    clustyp;                // Types of the association clusters (true - Type I; false - Type II)

		list<int>::iterator iti, itic;
		list<int>::iterator itass;
		list<bool>::iterator itasst;

		int  afterI = 0;
		for(iti = _currfac.anchors.begin(), itic = _currfac.cardinals.begin(); iti != _currfac.anchors.end(); iti++, itic++)
		{
			// Gathering all type II clusters before current type I
			int typeIoccur = *iti;
			if(typeIoccur > afterI)
			{
				assoc_clust.push_back(afterI);
				clustyp.push_back( false );

				for(int i = afterI+1; i < typeIoccur; i++)
				{
					if( this->weight[i] > this->weight[i-1] )
					{
						assoc_clust.push_back(i);
						clustyp.push_back( false );
					}
				}

			}

			// Current type I cluster
			assoc_clust.push_back(typeIoccur);
			clustyp.push_back( true );

			afterI = typeIoccur + *itic;
		}

		// Gathering all type II clusters after the last type I
		if(afterI < this->num)
		{
			assoc_clust.push_back(afterI);
			clustyp.push_back( false );
		}
		for(int i = afterI + 1; i < this->num; i++)
		{
			if( this->weight[i] > this->weight[i-1] )
			{
				assoc_clust.push_back(i);
				clustyp.push_back( false );
			}
		}

		// The second step - Cartesian multiplication

		list<DataVector> found_vecs;   // ${\cal V}_{F*}$

		// Forming "associated sets" $X^j_{assoc}$
		list<list<int> > ass_sets;
		for(itass = assoc_clust.begin(), itasst = clustyp.begin(); itass != assoc_clust.end(); itass++, itasst++)
		{
			list<int> X_j_assoc; 
			
			if(*itasst)   // If type I
			{
				X_j_assoc.push_back(*itass);
			}
			else          // If type II
			{
				list<int>::iterator itassn = itass;
				itassn++; 

				if(itassn != assoc_clust.end())
				{
					for(int jj = *itass; jj < *itassn; jj++)
					{
						X_j_assoc.push_back(jj);
					}
				}
				else
				{
					for(int jj = *itass; jj < this->num; jj++)
					{
						X_j_assoc.push_back(jj);
					}
				}
			}

			ass_sets.push_back(X_j_assoc); 
		}

		list<list<int> >::iterator itx;
		for(itx = ass_sets.begin(), itasst = clustyp.begin(); itx != ass_sets.end(); itx++, itasst++)
		{
			list<list<int> >::iterator itxn    = itx;
			list<bool>::iterator      itasstn = itasst;
			itxn++, itasstn++;

			// *itx = $X^j_assoc$, *itxn = $X^{j+1}_assoc$, *itasst and *itasstn - corresponding cluster types

			// Cartesian multiplication followed by the union
			if(itxn != ass_sets.end()/* && !(*itasst && *itasstn)*/)
			{
				// Cartesian
				list<Vector> part_found;

				list<int> X_j_assoc  = *itx;
				list<int> X_j1_assoc = *itxn;

				for(iti = X_j_assoc.begin(); iti != X_j_assoc.end(); iti++)
				{
					for(itic = X_j1_assoc.begin(); itic != X_j1_assoc.end(); itic++)
					{
						DataVector crit_vec(this->dim, *iti, *itic);

						crit_vec.start_typeI = *itasst;
						crit_vec.end_typeI   = *itasstn;

						crit_vec.FromPointToPoint(this->perm.points[_currfac.index_perm[*iti]], this->perm.points[_currfac.index_perm[*itic]]);
						
						// Union
						found_vecs.push_back(crit_vec);
					}
				}
			}
		}

		return found_vecs; 

	}


	list<Vector> GetCriticalVectors2(Facet _currfac)
	{
		list<Vector> crit_set;
		//
		//// Defining critical vectors connected with defining sets
		//for(int i = 0; i < _currfac.anchors.size(); i++)
		//{
		//	Vector cv(this->dim);
		//	// To the left
		//	if( _currfac.anchors[i]==1 || this->weight[_currfac.anchors[i]-1] > this->weight[_currfac.anchors[i]-2] )
		//	{
		//		cv.comb[0] = anchors[i]-1;
		//		cv.comb[1] = anchors[i];

		//		cv.FromPointToPoint(this->perm[indexperm[_currfac.anchors[i]]], this->perm[indexperm[_currfac.anchors[i]-1]]);

		//		crit_set.push_back(cv);
		//	}
		//	else
		//	{
		//		// Create vector for each point from the homogenous set to the left of curr_def_set

		//	}

		//	// To the right - the same
		//}
		//
		//// Defining critical vectors connected with emerging of new defining sets
		//// 1) Vectors between neighbors
		//// 2) Combintations from two neighboring "homogenous" sets (# = card(set_1) * card(set_2))


	}

protected:

	// To a given k-dimensional face find some orthogonal vector that is also orthogonal to its normal
	Vector FindRotDirection(kFace _facet)
	{
		// Defining some rotation vector 
		vector<vector<float> > def_vectors = GetDefVectors_kFace(_facet);

		// Augmenting the matrix:
		// 1. Orthogonal to normvec
		def_vectors[_facet.K()] = _facet.normalvec.coord;
		//def_vectors[j] = normvec.coord;
		// 2. Last coordinates = 0 except of the last = 1
		for(int k = _facet.K()+1; k < dim; k++)
		{
			vector<float> onevec(this->dim, 0);
			onevec[k] = 1;

			def_vectors[k] = onevec;
		}

		Vector rotvec(this->dim);

		vector<float> bf(this->dim, 0);
		bf[this->dim-1] = 1;
		rotvec.coord = _facet.SolveLinearSystem(def_vectors, bf);

		rotvec.Normalize();

		return rotvec;

	}

	// Find the neighboring (k+1)-face in a given direction
	kFace RemoveOneDoF(kFace _currface, Vector _rotvec)
	{
		kFace newface  = _currface;
		Vector normvec = _currface.normalvec;

		list<DataVector> crit_vectors = this->GetCriticalVectors(newface);

		list<DataVector>::iterator itcrit;
		float minplambda = 1000001;
		float minnlambda = 1;
		Vector deltavec(this->dim);

		DataVector minp, minn;

		// Searching for a minimal rotation on an angle < 90 degrees, if there is
		for(itcrit = crit_vectors.begin(); itcrit != crit_vectors.end(); itcrit++)
		{
			deltavec.coord = itcrit->coord;
			
			float lambda = - Vector::ScalarMultiply(deltavec, normvec) / Vector::ScalarMultiply(deltavec, _rotvec);

			if(lambda > 0.0000001f)
			{
				if (lambda < minplambda)
				{
					minplambda = lambda;
					minp       = *itcrit;
				}								
			}
		}

		// If the rotation should done on the angle > 90 degrees
		if( !minp.defined() )
		{
			for(itcrit = crit_vectors.begin(); itcrit != crit_vectors.end(); itcrit++)
			{
				deltavec.coord = itcrit->coord;
				
				float lambda = - Vector::ScalarMultiply(deltavec, normvec) / Vector::ScalarMultiply(deltavec, _rotvec);
				
				if(lambda <= -0.0000001f)
				{
					if (lambda < minnlambda)
					{
						minnlambda = lambda;
						minn       = *itcrit;
					}
				}
			}
		}

		DataVector neighn;      // Found neighbor

		if( minp.defined() )
		{
			neighn  = minp;

			Vector vecp = normvec;
			Vector deltap = _rotvec;
			deltap.Scale(minplambda);
			vecp.Add(deltap);

			newface.normalvec = vecp;
		}
		else
		{
			neighn = minn;

			Vector vecn = normvec;
			vecn.Reverse();
			Vector deltan = _rotvec;
			deltan.Scale(-minnlambda);
			vecn.Add(deltan);

			newface.normalvec = vecn;
		}


		// Normalizing the normalvec (normvec & _rotvec must have length=1)
		float tonormal;
		if( minp.defined() )
			tonormal = 1.0f / (sqrt(1.0f + minplambda * minplambda));
		else
			tonormal = 1.0f / (sqrt(1.0f + minnlambda * minnlambda));

		newface.normalvec.Scale(tonormal);

		// If 2 Type II clusters do associate, create a new ${\cal{A}}_q$ 

		// Modifying ${\cal X}_F$ by adding the found vector (to some existing or new ${\cal{A}}_q$)

		// Editing index_perm
		// Warning: let neighn - the found vector to be intoduced

		if(neighn.start > neighn.end)
		{
			swap(neighn.start, neighn.end);

			swap(neighn.start_typeI, neighn.end_typeI);
		}

		list<int>::iterator srcit1, srcit2; 
		
		// Type I - Type I
		if(neighn.start_typeI && neighn.end_typeI)
		{
			// newface.index_perm stays unchanged

			// Locating corresponding ${\cal{A}}$s
			for(srcit1 = newface.anchors.begin(), srcit2 = newface.cardinals.begin();
				srcit1 != newface.anchors.end(); 
				srcit1++, srcit2++)
			{
				if(*srcit1 > neighn.start) break;
			}
			
			newface.anchors.erase(srcit1); 
			
			int added = *srcit2;
			srcit2--;
			*srcit2 += added;
			srcit2++;
			newface.cardinals.erase(srcit2);
			
		}
		// Type II - Type II
		else if(!neighn.start_typeI && !neighn.end_typeI)
		{
			// Shifting start point to the right border of the type II cluster
			int qq = neighn.start;
			while(this->weight[qq+1] == this->weight[qq])
			{
				qq++;
			}

			if(qq > neighn.start)
			{
				swap(newface.index_perm[qq], newface.index_perm[neighn.start]);
				neighn.start = qq;
			}

			// Shifting end point to the left border of the type II cluster
			qq = neighn.end;
			while(this->weight[qq-1] == this->weight[qq])
			{
				qq--;
			}

			if(qq < neighn.end)
			{
				swap(newface.index_perm[qq], newface.index_perm[neighn.end]);
				neighn.end = qq;
			}

			// newface.index_perm is ready

			// Locating place for a new ${\cal{A}}$
			for(srcit1 = newface.anchors.begin(), srcit2 = newface.cardinals.begin();
				srcit1 != newface.anchors.end(); 
				srcit1++, srcit2++)
			{
				if(*srcit1 > neighn.start) break;
			}

			newface.anchors.insert(srcit1, neighn.start);
			newface.cardinals.insert(srcit2, 2);
		}
		else
		{
			// Type I - Type II
			if(neighn.start_typeI)
			{
				//// Shifting end point to the left border of the type II cluster
				//int qq = neighn.end;
				//while(this->weight[qq-1] == this->weight[qq])
				//{
				//	qq--;
				//}

				// newface.index_perm is ready

				// Locating the ${\cal{A}}$
				for(srcit1 = newface.anchors.begin(), srcit2 = newface.cardinals.begin();
					srcit1 != newface.anchors.end(); 
					srcit1++, srcit2++)
				{
					if(*srcit1 > neighn.start) break;
				}

				srcit2--;
				*srcit2 += 1;
				srcit1--;

				if(*srcit1 + *srcit2 - 1 < neighn.end)
				{
					swap(newface.index_perm[*srcit1 + *srcit2 - 1], newface.index_perm[neighn.end]);
					neighn.end = *srcit1 + *srcit2;
				}
			}
			// Type II - Type I
			else
			{
				//// Shifting start point to the right border of the type II cluster
				//int qq = neighn.start;
				//while(this->weight[qq+1] == this->weight[qq])
				//{
				//	qq++;
				//}

				// newface.index_perm is ready

				// Locating the ${\cal{A}}$
				for(srcit1 = newface.anchors.begin(), srcit2 = newface.cardinals.begin();
					srcit1 != newface.anchors.end(); 
					srcit1++, srcit2++)
				{
					if(*srcit1 > neighn.start) break;
				}

				*srcit1 -= 1;
				*srcit2 += 1;

				if(*srcit1 > neighn.start)
				{
					swap(newface.index_perm[*srcit1], newface.index_perm[neighn.start]);
					neighn.start = *srcit1;
				}
			}
		}

		newface.CatchOneDoF();
		return newface;
	}

public:

	float IntersectionCriterion(Point _orig, Point _mu, Facet _fac)
	{
		float crit = -_fac.abs_member / ( Vector::ScalarMultiply(_fac.normalvec, Vector::CreateVector(_mu)) ); 

		return crit;
	}

	// The method for searching the facet being intersected by the line passing through the point _orig (default - the origin)
	// and with direction given by the vector _phi
        // Special cases: // 0 - normal; 1 - unlimited; 2 - no solution
	Facet FindIntersectedFacet(Point _orig, Point _mu, int& _case)
	{
		Vector phi(this->dim);    // The line between _orig and _mu
		phi.FromPointToPoint(_mu,_orig);

		// ************************* Smart first facet ****************************


		// Create an "empty" first facet
		kFace ffac(this->dim, 0);

		// Initial $\pi(i)$ - trivial
		vector<int>  initial_ind(num);
		for(int fu = 0; fu < num; fu++)
			initial_ind[fu] = fu;

		ffac.index_perm = initial_ind;

		// Getting the first permutation - only once for the first facet!
		this->perm.Support(phi);

		// Now, only anchors and cardinals must be constructed (also normalvec)

		ffac.normalvec = phi;
		ffac.normalvec.Normalize();

		// Grasping d-1 degrees of freedom
		for(int j = 0; j < this->dim - 1; j++)
		{
			// Create the direction to shift the normal vector
			Vector rotvec = this->FindRotDirection(ffac);

			// Moving normvec by means of rotvec, which is orthogonal to normvec

			ffac = this->RemoveOneDoF(ffac, rotvec);
	
		}
		
		//ffac.normalvec.Normalize();
		Facet found_facet(ffac);
		found_facet.doubled   = false;

		Point ffnode(this->dim);
		ffnode.coord = CurrentExtremePoint(found_facet.index_perm);
		found_facet.CalculateAbsoluteMemberH(ffnode);
		
		// Take found_facet as a start point for searching the intersection
		Facet currfacet = found_facet;

		this->trimmed_region.push_back(currfacet);


		::RecordTime();
		
		// ************************************************************************

		// *********** Start travelling towards the intersection ******************
				
		bool intersection_found;
                bool front_achieved = false;
                float curr_crit_back = 5555555, curr_crit = -55;

                float init_crit = this->IntersectionCriterion(_orig, _mu, currfacet);

                if(currfacet.abs_member >= 0)
                {
                        curr_crit   = init_crit;
                        front_achieved = true; 
                }
                else if(init_crit >= 0)
                    curr_crit_back = init_crit;

                do
		{
			intersection_found = true;
                        float currfac_crit = this->IntersectionCriterion(_orig, _mu, currfacet);

			Facet nicer_facet = currfacet; 

			// Removing one vector thus enabling rotating

			// Iterating through ${\cal{A}}_l \in {\cal X}_F$
			list<int>::iterator iterset, itercard;
			for(iterset = currfacet.anchors.begin(), itercard = currfacet.cardinals.begin(); 
				iterset != currfacet.anchors.end(); 
				iterset++, itercard++)
			{
				// "Cracking" each ${\cal{A}}_l$ into 2 sets. There are $2^{\|{\cal{A}}\|}-2$

				list<boost::dynamic_bitset<> > crackmasks;

				for(int ccc = 1; ccc < pow(2.0, *itercard) - 1.5; ccc++ )
				{
					boost::dynamic_bitset<> crackmask(*itercard, ccc);

					crackmasks.push_back(crackmask);
				}

				// Filtering crackmaps according to the condition (3)

				list<boost::dynamic_bitset<> >::iterator itcrack; 

				//for(int k = *iterset; k < *iterset+*itercard; k++)
				for(itcrack = crackmasks.begin(); itcrack != crackmasks.end(); itcrack++)
				{
					int bulk = itcrack->count();

					// Moving in the permutation according to the cracking map
					if( (this->weight[*iterset+*itercard-1] > this->weight[*iterset+bulk]) && ( this->weight[*iterset+bulk-1] > this->weight[*iterset]) ||
						(bulk == 1)           && (this->weight[*iterset+*itercard-1] > this->weight[*iterset+bulk]) ||
						(bulk == *itercard-1) && (this->weight[*iterset+bulk-1]      > this->weight[*iterset]) ||
						*itercard == 2)
					{
						// Splitting the ${\cal{A}}_l$ according to the cracking map

						//::ResumeCumulTime();
						// Be cautious: currfacet must be recovered
						currfacet.anchors.insert(iterset, *iterset);
						int origset = *iterset;
						*iterset  = *iterset + bulk;
						list<int>::iterator itsetnew  = iterset;
						iterset--;

						currfacet.cardinals.insert(itercard, bulk);
						int origcard = *itercard;
						*itercard = *itercard - bulk;
						list<int>::iterator itcardnew = itercard;
						itercard--;

						// When removing only 1 point, or eliminating ${\cal{A}}_l$ containing only 2 points
						if(bulk == 1)
						{
							list<int>::iterator empset = iterset;
							list<int>::iterator empcard = itercard;
							iterset++; iterset++;      // 2-step move
							itercard++; itercard++;      // 2-step move
							currfacet.anchors.erase(empset);
							currfacet.cardinals.erase(empcard);
						}

						if(bulk == itcrack->size()-1)
						{
							currfacet.anchors.erase(itsetnew);
							currfacet.cardinals.erase(itcardnew);
						}
						
						kFace newridge = currfacet;
						newridge.RelaxOneDoF();
						
						// Modifying index_perm
						int ww1 = origset, ww2 = origset + origcard - 1;
						for(int ww = 0; ww < itcrack->size(); ww++)
						{
							if((*itcrack)[ww])
							{
								newridge.index_perm[ww1] = currfacet.index_perm[origset + ww];
								ww1++;
							}
							else
							{
								newridge.index_perm[ww2] = currfacet.index_perm[origset + ww];
								ww2--;
							}
						}

						// Taking a vector to detect direction of an infinitesimal rotation (either clockwise or counterclockwise)
						Vector infinitesimalvec(this->dim);
						infinitesimalvec.FromPointToPoint(perm.points[newridge.index_perm[origset]], perm.points[newridge.index_perm[origset + origcard - 1]]);

						// Saving removed vector (one of the possible)
						//newridge.old_vector.FromPointToPoint(this->perm.points[newridge.index_perm[origset]], this->perm.points[newridge.index_perm[origset + bulk]]);
						
						// Recovering currfacet
						if(bulk < itcrack->size()-1)
						{
							currfacet.anchors.erase(itsetnew);
							currfacet.cardinals.erase(itcardnew);
						}

						if(bulk == 1)
						{
							currfacet.anchors.insert(iterset, origset);
							currfacet.cardinals.insert(itercard, origcard);
							iterset--;
							itercard--;
						}

						*iterset  = origset;
						*itercard = origcard;

						// currfacet recovered

						Vector rotvec = this->FindRotDirection(newridge);

						if( Vector::ScalarMultiply(infinitesimalvec, rotvec) <  -0.00001f)
							rotvec.Reverse();

						// Now rotvec is orthogonal to currfacet.normalvec. <rotvec;currfacet.normalvec> is a B2 basis
						// Direction of rotation (clcws. or counterclcws.) is shown by infinitesimalvec 

						Facet neigh_facet(this->RemoveOneDoF(newridge, rotvec));
						neigh_facet.marked = false;
						Point idnode(this->dim);
						idnode.coord = CurrentExtremePoint(neigh_facet.index_perm);
						neigh_facet.CalculateAbsoluteMemberH(idnode);

						// Maximizing the criterion
						float neigh_crit = this->IntersectionCriterion(_orig, _mu, neigh_facet);
                                                
						// Note that we are interested in the intersections only on the lower boundary,
						// that is under the expectation.
                                                // The intersection with the back is interesting in detecting special cases:
                                                // 1. Only back - no solution. 2. Back earlier - unlimited solution. 

                                                if(currfacet.abs_member >= 0)
                                                {
                                                    // Front - back
                                                    if(neigh_facet.abs_member < 0)
                                                    {
                                                        if(neigh_crit < curr_crit_back && neigh_crit >= 0 )
                                                        {
                                                            curr_crit_back   = neigh_crit;
                                                        }
                                                    }
                                                    // Front - front
                                                    else
                                                    {
                                                        if(neigh_crit > curr_crit)
                                                        {
                                                            intersection_found = false;

                                                            nicer_facet = neigh_facet;
                                                            curr_crit   = neigh_crit;
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    // Back - back
                                                    if(neigh_facet.abs_member < 0)
                                                    {
                                                        if(!front_achieved && fabs(neigh_crit) < fabs(currfac_crit))
                                                        {
                                                            intersection_found = false;
                                                            nicer_facet = neigh_facet;
                                                        }
                                                    }
                                                    // Back - front
                                                    else
                                                    {
                                                        intersection_found = false;

                                                        // Just set as currfacet and break
                                                        front_achieved = true;
                                                        nicer_facet = neigh_facet;
                                                    }
                                                }

					} // end processing a feasible crackmask
				}  // end of loop for a current crackmask

			}

                        if(curr_crit_back < curr_crit)
                        {
                            _case = 1;   // Unlimited solution
                            break;
                        }

			currfacet = nicer_facet;

			this->trimmed_region.push_back(currfacet);

		}while(!intersection_found);

                if(curr_crit < 0) _case = 2; // No solution - the efficient set contains the origin

                return currfacet;
                
	}


private:

	// Returns a set of splittings of the _sublist into a selected combination and the rest 
	template<class ArbClass> list<vector<vector<ArbClass> > > SelectCombins(vector<ArbClass> _vec, int _num)
	{
		int veclen = _vec.size(); 
		int numofcombs = Comb(veclen, _num);

		list<vector<vector<ArbClass> > > reslist;

		vector<int> comb(_num), contracomb(veclen - _num);
		for(int i=0; i < _num; i++)
		{
			comb[i]       = i;
		}
		for(int i=_num; i < veclen; i++)
		{
			contracomb[i-_num] = i;
		}

		vector<vector<ArbClass> > splitt(2);
		vector<ArbClass> empt1(_num);          splitt[0] = empt1;
		vector<ArbClass> empt2(veclen - _num); splitt[1] = empt2;

		for(int w = 0; w < _num; w++)
			(splitt[0])[w] = _vec[comb[w]]; 
		for(int w = 0; w < veclen - _num; w++)
			(splitt[1])[w] = _vec[contracomb[w]];
		//splitt[0] = comb; splitt[1] = contracomb;

		reslist.push_back(splitt);
		
		for(int q=0; q < numofcombs-1; q++)
		{
			comb = NextCombination(comb, veclen);
			
			int hund = 0;
			for(int i=0; i < veclen; i++)
			{
				if(hund < _num && comb[hund] == i )
					hund++;
				else contracomb[i-hund] = i; 
			} 
		
			for(int w = 0; w < _num; w++)
				(splitt[0])[w] = _vec[comb[w]]; 
			for(int w = 0; w < veclen - _num; w++)
				(splitt[1])[w] = _vec[contracomb[w]];
			//splitt[0] = comb; splitt[1] = contracomb;
		
			reslist.push_back(splitt);
		}

		return reslist;
	}

public:
	
	void CalculateAllVertices()
	{
		list<Facet>::iterator facetit;

		// For each facet
		for(facetit = this->trimmed_region.begin(); facetit != this->trimmed_region.end(); facetit++)
		{
			list<int>::iterator iterit, itercard;

			// Combinations within ${\cal A}_l$ are independent of ones in ${\cal A}_k, k \ne l$
			vector<list<vector<vector<int> > > > als;
			for(iterit = facetit->anchors.begin(), itercard = facetit->cardinals.begin();
				iterit != facetit->anchors.end();
				iterit++, itercard++)
			{
				// Create a list from $A_l$ and weights homogeneity map
				vector<int> al, alordered;
				list<int> alwmap;
				int stopal = *iterit+*itercard;
					
				int lastjump = *iterit;
				for(int i = *iterit; i < stopal; i++)
				{
					//al.push_back(facetit->index_perm[i]);
						
					if( (i > *iterit) && (weight[i] > weight[i-1]) )
					{
						alwmap.push_back(i - lastjump);
						lastjump = i; 
					}
				}
				alwmap.push_back(stopal - lastjump);

				al.assign(facetit->index_perm.begin() + *iterit, facetit->index_perm.begin() + stopal);

				list<vector<vector<int> > > alcombs; // List of 2-element vectors - splittings of A_l into permuted and unpermuted parts
				vector<vector<int> > defcomb(2);
				defcomb[0] = alordered; defcomb[1] = al;
				alcombs.push_back(defcomb);

				list<vector<vector<int> > >::iterator itspl,itpartspl;
				do
				{
					list<vector<vector<int> > > nextalcombs;

					for(itspl = alcombs.begin(); itspl != alcombs.end(); itspl++)
					{
						list<vector<vector<int> > > partcombs = SelectCombins((*itspl)[1], alwmap.front());
						
						// Augmenting partial combination with an already ordered part
						for(itpartspl = partcombs.begin(); itpartspl != partcombs.end(); itpartspl++)
						{
							(*itpartspl)[0].insert((*itpartspl)[0].begin(), (*itspl)[0].begin(), (*itspl)[0].end());
						}

						nextalcombs.splice(nextalcombs.begin(), partcombs);
					}

					alcombs = nextalcombs;

					alwmap.pop_front();
				}while(!alwmap.empty());

				als.push_back(alcombs);
			}

			// Generating final permutations (combining independent variants from all $A_l$)
			list<vector<int> > ready_perms;
			vector<list<vector<vector<int> > > >::iterator iterals;

			ready_perms.push_back(facetit->index_perm);

			for(iterit = facetit->anchors.begin(), itercard = facetit->cardinals.begin(), iterals = als.begin();
				iterit != facetit->anchors.end();
				iterit++, itercard++, iterals++)
			{
				list<vector<int> > step_ready_perms;
				int stopal = *iterit+*itercard;

				for(list<vector<int> >::iterator  itsub = ready_perms.begin(); itsub != ready_perms.end(); itsub++)
				{
					for(list<vector<vector<int> > >::iterator tmpal = iterals->begin(); tmpal != iterals->end(); tmpal++)
					{
						vector<int> augmperm = *itsub;
						
						for(int j = *iterit; j < stopal; j++)
						{
							augmperm[j] = ((*tmpal)[0])[j-*iterit];
						}

						step_ready_perms.push_back(augmperm);
					}
				}

				ready_perms = step_ready_perms;
			}

			// Final genearing of vertices

			facetit->nodes.clear();

			for(list<vector<int> >::iterator  itsub = ready_perms.begin(); itsub != ready_perms.end(); itsub++)
			{
				ExtremePoint newnode(this->dim);

				newnode.coord = CurrentExtremePoint(*itsub);
				facetit->nodes.push_back(newnode);
			}
			
		}
	}
	
	void CalculateAllVertices3d()
	{
		if(dim == 3)
		{
			list<Facet>::iterator facetit;

			// For each facet
			for(facetit = this->trimmed_region.begin(); facetit != this->trimmed_region.end(); facetit++)
			{
				// If there are 2 Al's
				if(facetit->anchors.size() == 2)
				{
					// Faces with 4 vertices are generated
					facetit->truncated = true;

					vector<int> piperm = facetit->index_perm;

					ExtremePoint idnode(this->dim);

					idnode.coord = CurrentExtremePoint(piperm);
					facetit->nodes.push_back(idnode);

					swap(piperm[facetit->anchors.front()], piperm[facetit->anchors.front()+1]);
					idnode.coord = CurrentExtremePoint(piperm);
					facetit->nodes.push_back(idnode);

					swap(piperm[facetit->anchors.back()], piperm[facetit->anchors.back()+1]);
					idnode.coord = CurrentExtremePoint(piperm);
					facetit->nodes.push_back(idnode);

					swap(piperm[facetit->anchors.front()], piperm[facetit->anchors.front()+1]);
					idnode.coord = CurrentExtremePoint(piperm);
					facetit->nodes.push_back(idnode);

				}
				// If there is only 1 Al
				else
				{
					int unianc = facetit->anchors.front();

					// If not all weights are different - 3 vertices
					if( this->weight[unianc] == this->weight[unianc + 1])
					{
						facetit->truncated = false;

						vector<int> piperm = facetit->index_perm;

						ExtremePoint idnode(this->dim);

						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
						
						swap(piperm[unianc], piperm[unianc + 2]);
						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
						
						swap(piperm[unianc + 1], piperm[unianc + 2]);
						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
					}
					else if( this->weight[unianc+1] == this->weight[unianc + 2] )
					{
						facetit->truncated = false;

						vector<int> piperm = facetit->index_perm;

						ExtremePoint idnode(this->dim);

						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
						
						swap(piperm[unianc], piperm[unianc + 2]);
						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
						
						swap(piperm[unianc], piperm[unianc + 1]);
						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
					}
					// Else - 6 vertices
					else
					{
						facetit->truncated = true;

						vector<int> piperm = facetit->index_perm;

						ExtremePoint idnode(this->dim);

						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
						
						swap(piperm[unianc], piperm[unianc + 1]);
						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
						
						swap(piperm[unianc + 1], piperm[unianc + 2]);
						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
						
						swap(piperm[unianc], piperm[unianc + 1]);
						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
						
						swap(piperm[unianc + 1], piperm[unianc + 2]);
						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
						
						swap(piperm[unianc], piperm[unianc + 1]);
						idnode.coord = CurrentExtremePoint(piperm);
						facetit->nodes.push_back(idnode);
					}
					
					/*list<ExtremePoint>::iterator iterex;
					for(iterex = facetit->nodes.begin(); iterex != facetit->nodes.end(); iterex++)
					{
						Vector vec(this->dim);
						vec.coord = iterex->coord;
						float chdelta = fabs(Vector::ScalarMultiply(vec, facetit->normalvec) + facetit->abs_member);
						if(chdelta > 0.0001f)
						{
							vec.Reverse(); // ERROR
						}
					}*/
				}
			}
		}

	}

	list<Facet> ReadyTR()
	{
		return this->trimmed_region;
	}
};


// ******************************* Class for receiving original data cloud **************************

class InputAg : public Agent
{
public:

	int   dim  ;                  // Dimension of the cloud
	vector<Point> cloud;
	Point centroid;

	vector<float> cvec;          // c-vector - the objective
        float bval;                  // b - the value of the right-hand side

private:

	string WMTD_type;
	float depth;                  // Parameter for computing trimmed region
	int   num  ;                  // Number of points
	vector<float> manweights;

public:

	InputAg():
	  manweights(),
                  cvec()
	{
		dim   = 3;
		num   = 0;

		depth = 1;
	}

	// Receive cloud of points from the specified text file
	int Receive(char* _source, char* _dir)
	{
		try
		{
			ifstream is(strcat(_dir,_source));     // Open file with data

			is >> dim;

			for(int j=0; j<dim; j++)
			{
				cvec.push_back(j);
				is >> cvec[j];
			}

                        is >> bval;

			is >> WMTD_type;

			is >> depth;

			is >> num;

			// Reading coordinates to objects...
			for(int i=0; i<num; i++)
			{
				Point p(dim);
				for(int j=0; j<dim; j++)
				{
					is >> p.coord[j];
				}
				
				cloud.push_back(p);       // Adding a point to cloud

			}
			
			is.close();

			return 0;
		}
		catch(std::exception ex)
		{
			return -1;
		}
		
		return 0;
	}

	vector<float> ProcessData(ResultAg* _result_agent)
	{
		// Creating main processor for realizing algorithm 
		ProcessAg* processor = new ProcessAg(WMTD_type, depth, dim, num, cloud, manweights);

		this->cloud    = processor->perm.points;
		this->centroid = processor->initial_center;

		// Default origin - the point \bm0
		Point orignull(this->dim);
		vector<float> nullvec(dim, 0);
		orignull.coord = nullvec;

		::StartTiming();
		//return processor->Compute();
		Point cray = orignull; cray.coord = this->cvec;
                int special_case = 0; // 0 - normal; 1 - unlimited; 2 - no solution
		Facet intersected_facet = processor->FindIntersectedFacet(orignull, cray, special_case);
		::RecordTime();

                // 0 - normal; 1 - unlimited; 2 - no solution
                if(special_case == 1)
                {
                    vector<float> answ(intersected_facet.normalvec.coord.size(),-11);
                    return answ;
                }
                else if(special_case == 2)
                {
                    vector<float> answ(intersected_facet.normalvec.coord.size(),-22);
                    return answ;
                }

		// Normalizing the inverted vector to the component sum of 1
		vector<float> optplio= intersected_facet.normalvec.coord;

                float scaling = this->bval / intersected_facet.abs_member;
		for(int i = 0; i < optplio.size(); i++)
		{
			optplio[i] = -optplio[i] * scaling;
		}
		/*
                float sum_of_components = 0;
		for(int i = 0; i < optplio.size(); i++)
		{
			//optplio[i] = 1.0f / optplio[i];
			sum_of_components += optplio[i];
		}
		for(int i = 0; i < optplio.size(); i++)
		{
			optplio[i] = optplio[i] / sum_of_components;
		}
                 */


		if(this->dim == 3)
			processor->CalculateAllVertices3d();
		/*else
			processor->CalculateAllVertices();*/

		//return processor->ReadyTR();
		_result_agent->trimmed_region = processor->ReadyTR();

		return optplio;

	}

};

// ************************************************************************************************
// *********************** Main Function **********************************************************

extern "C"{

int  SolveSLP_WMTR(char** nameofsource, char** _wdir)			// Window Show State
{
	
	bool single_mode = true;
	
	ResultAg* result_agent = new ResultAg();
	InputAg*  input_agent  = new InputAg();
	
	int res = input_agent->Receive(*nameofsource, *_wdir);

	vector<float> optportfac     = input_agent->ProcessData(result_agent);
	result_agent->data_cloud     = input_agent->cloud;
	result_agent->coord_center   = input_agent->centroid;
	result_agent->coord_center.coord = input_agent->cvec;

	// Print results
	result_agent->PrintSolution(optportfac);
	result_agent->GLSceneToExport();
}
}
