#include "network.h"
#include "random.h"

void Network::resize(const size_t &n, double inhib) {
    size_t old = size();
    neurons.resize(n);
    if (n <= old) return;
    size_t nfs(inhib*(n-old)+.5);
    set_default_params({{"FS", nfs}}, old);
}

void Network::set_default_params(const std::map<std::string, size_t> &types,
                                 const size_t start) {
    size_t k(0), ssize(size()-start), kmax(0);
    std::vector<double> noise(ssize);
    _RNG->uniform_double(noise);
    for (auto I : types) 
        if (Neuron::type_exists(I.first)) 
            for (kmax+=I.second; k<kmax && k<ssize; k++) 
                neurons[start+k].set_default_params(I.first, noise[k]);
    for (; k<ssize; k++) neurons[start+k].set_default_params("RS", noise[k]);
}

void Network::set_types_params(const std::vector<std::string> &_types,
                               const std::vector<NeuronParams> &_par,
                               const size_t start) {
    for (size_t k=0; k<_par.size(); k++) {
        neurons[start+k].set_type(_types[k]);
        neurons[start+k].set_params(_par[k]);
    }
}

void Network::set_values(const std::vector<double> &_poten, const size_t start) {
    for (size_t k=0; k<_poten.size(); k++) 
        neurons[start+k].potential(_poten[k]);
}

bool Network::add_link(const size_t &a, const size_t &b, double str) {
    if (a==b || a>=size() || b>=size() || str<1e-6) return false;
    if (links.count({a,b})) return false;
    if (neurons[b].is_inhibitory()) str *= -2.0;
    links.insert({{a,b}, str});
    return true;
}

size_t Network::random_connect(const double &mean_deg, const double &mean_streng) {
    links.clear();
    std::vector<int> degrees(size());
    _RNG->poisson(degrees, mean_deg);
    size_t num_links = 0;
    std::vector<size_t> nodeidx(size());
    std::iota(nodeidx.begin(), nodeidx.end(), 0);
    for (size_t node=0; node<size(); node++) {
        _RNG->shuffle(nodeidx);
        std::vector<double> strength(degrees[node]);
        _RNG->uniform_double(strength, 1e-6, 2*mean_streng);
        int nl = 0;
        for (size_t nn=0; nn<size() && nl<degrees[node]; nn++)
            if (add_link(node, nodeidx[nn], strength[nl])) nl++;
        num_links += nl;
    }
    return num_links;
}

std::vector<double> Network::potentials() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].potential());
    return vals;
}

std::vector<double> Network::recoveries() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].recovery());
    return vals;
}

void Network::print_params(std::ostream *_out) {
    (*_out) << "Type\ta\tb\tc\td\tInhibitory\tdegree\tvalence" << std::endl;
    for (size_t nn=0; nn<size(); nn++) {
        std::pair<size_t, double> dI = degree(nn);
        (*_out) << neurons[nn].formatted_params() 
                << '\t' << dI.first << '\t' << dI.second
                << std::endl;
    }
}

void Network::print_head(const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons)
            if (In.is_type(It.first)) {
                (*_out) << '\t' << It.first << ".v"
                        << '\t' << It.first << ".u"
                        << '\t' << It.first << ".I";
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << "RS.v" << '\t' << "RS.u" << '\t' << "RS.I";
                break;
            }
    (*_out) << std::endl;
}

void Network::print_traj(const int time, const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    (*_out)  << time;
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons) 
            if (In.is_type(It.first)) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    (*_out) << std::endl;
}




////////////////////////////////////////////////////////


/*
 Finds the list of neurons with incoming connections to \p n.
  \param n : the index of the receiving neuron.
  \return a vector of pairs {neuron index, link intensity}.
*/
std::vector<std::pair<size_t, double> > Network::neighbors(const size_t& n) const
{
	///on regarde dans tableau links tous les neurones connectés à n
	///on les ajoute à vector qu'on retourne
	
	///typedef std::map<std::pair<size_t, size_t>, double> linkmap;
	///map<pair<neurone1, neurone2>, intensité>, neurone1 reçoit intensité de neurone2
	std::vector<std::pair<size_t, double> > voisins;
	std::pair<size_t, double> voisin;
	
	///links.begin() pas bonne option parce qu'on parcourt toute la liste de neurones --> trop long, mieux links.lower_bound()
	///cend au lieu de end quand fn const
	///pour marquer la fin : l'indice change, on passe à autre neurone dans links
	for(auto i = links.lower_bound({n, 0}); i != links.cend() and (i->first).first == n; ++i)
	{
		voisin.first = (i->first).second;
		voisin.second = i->second;
		voisins.push_back(voisin);
	}
	
	return voisins;
}


/*
Calculates the number and total intensity of connections to neuron \p n.
  \param n : the index of the receiving neuron.
  \return a pair {number of connections, sum of link intensities}.
*/

std::pair<size_t, double> Network::degree(const size_t& n) const
{
	std::pair<size_t, double> intensite;
	double sommeIntensites(0.0);
	
	///on doit appeler fn neighbors pour voir neurones connectés, prendre leur intensité, additionner
	///tableau de size_t,double
	for(size_t i=0 ; i<neighbors(n).size() ; ++i)
	{
		sommeIntensites += neighbors(n)[i].second;
	}
	
	intensite.first = neighbors(n).size();
	intensite.second = sommeIntensites;
	
	return intensite;
}


/*
Performs one time-step of the simulation.
  \param input : a vector of random values as thalamic input, one value for each neuron. 
  * The variance of these values corresponds to excitatory neurons.
  \return the indices of firing neurons.
*/

std::set<size_t> Network::step(const std::vector<double>& thalamicInput)
{
	///on prend neurones voisins, on regarde si firing ou pas, puis si inhibiteur ou excitateur
	
	std::set<size_t> firingNeurons;
	
	///on regarde chaque neurone
	for(size_t i=0; i<neurons.size(); ++i)
	{
		std::vector<std::pair<size_t, double> > voisins(neighbors(i));
		double sommeInhibiteur(0.0);
		double sommeExcitateur(0.0);
		
		///on regarde les voisins du neurone
		for(size_t j=0; j<voisins.size(); ++j)
		{
			///on regarde si voisin est firing
			if(neurons[voisins[j].first].firing())
			{
				///on regarde si voisin firing est inhibiteur
				if(neurons[voisins[j].first].is_inhibitory())
				{
					///on ajoute intensité à somme inhibiteur
					sommeInhibiteur += voisins[j].second;
				} else {
					///on ajoute intensité à somme excitateur
					sommeExcitateur += voisins[j].second;
				}
			}
		}
		
		if(neurons[i].is_inhibitory())
		{
			neurons[i].input(thalamicInput[i]*2.0/5.0 + 0.5*sommeExcitateur + sommeInhibiteur);
		} else {
			neurons[i].input(thalamicInput[i] + 0.5*sommeExcitateur + sommeInhibiteur);
		}
		
		if(neurons[i].firing())
		{
			firingNeurons.insert(i);
			neurons[i].reset();
		} else {
			neurons[i].step();
		}
	}
	
	return firingNeurons;
}






