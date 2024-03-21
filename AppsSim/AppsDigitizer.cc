#include "AppsDigitizer.hh"
#include "G4AutoLock.hh"
#include "G4Threading.hh"
#include <chrono>
#include "G4SystemOfUnits.hh"
#include "AppsInput.hh"
#include <unistd.h>
#include <tuple>
#include <atomic>
#include <mutex>

std::mutex countMutex;
std::atomic<long long> EventCount = 0; // variable that stores the current time of the simulation

using namespace std;

AppsDigitizer::AppsDigitizer()
{

	// 10000 is arbitrary for now
	StoredEvents.reserve(10000);
}

AppsDigitizer::~AppsDigitizer()
{
}

void AppsDigitizer::AddEvent(G4double time, G4double energy, G4int pileupflag)
{

	EventCount++;

	///	cout << " Event count: " << EventCount << " "  << StoredEvents.size() << " " <<  G4Threading::G4GetThreadId() << endl;

	if (EventCount > 1000)
	{
		//	auto start = std::chrono::high_resolution_clock::now();
		//		cout << " ------------- Build time ---------- for a number of:  " << EventCount << " size of " << StoredEvents.size() << " by " <<  G4Threading::G4GetThreadId() << endl;
		EventCount = 0;

		countMutex.lock();

		cout << endl;

		BuildDigitalEvent();

		countMutex.unlock();

		auto finish = std::chrono::high_resolution_clock::now();
		//	  std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count() << "ns\n";
	}

	// cout << StoredEvents.size() << endl;

	// cout << time/CLHEP::ns << " " << energy << " " << G4Threading::G4GetThreadId() << endl;

	StoredEvents.push_back(make_tuple(time, energy, pileupflag));
}

void AppsDigitizer::BuildDigitalEvent()
{

	//	cout << " in the build size " << StoredEvents.size() << endl;

	AppsInput *FInput = AppsInput::Instance();

	if (FInput->GetBeamType() == "GBS")
	{

		sort(StoredEvents.begin(), StoredEvents.end());

		cout << endl;
		cout << "Initial list " << endl;

		for (unsigned int i = 0; i < StoredEvents.size(); i++)
		{

			cout << get<0>(StoredEvents[i]) << " " << get<1>(StoredEvents[i]) << " " << get<2>(StoredEvents[i]) << " " << G4Threading::G4GetThreadId() << endl;
		}

		tbb::concurrent_vector<tuple<double, double, int>> RisePileUp;
		tbb::concurrent_vector<tuple<double, double, int>> TailPileUp;

		for (unsigned int i = 0; i < StoredEvents.size(); i++)
		{

			double currentTime = get<0>(StoredEvents[i]);
			double currentEnergy = get<1>(StoredEvents[i]);
			int pileUpFlag = get<2>(StoredEvents[i]);

			while (i < StoredEvents.size() - 1)
			{

				if ((get<0>(StoredEvents[i + 1]) - currentTime) < 200)
				{

					currentTime = (currentTime + get<0>(StoredEvents[i + 1])) / 2;
					currentEnergy += get<1>(StoredEvents[i + 1]);
					pileUpFlag += get<2>(StoredEvents[i + 1]);
				}
				else
				{

					break;
				}

				i++;
			}

			RisePileUp.push_back(make_tuple(currentTime, currentEnergy, pileUpFlag));
		}

		cout << endl;
		cout << "Intermediate list " << endl;

		for (unsigned int i = 0; i < RisePileUp.size(); i++)
		{

			cout << get<0>(RisePileUp[i]) << " " << get<1>(RisePileUp[i]) << " " << get<2>(RisePileUp[i]) << " " << G4Threading::G4GetThreadId() << endl;
		}

		for (unsigned int i = 0; i < RisePileUp.size(); i++)
		{

			double currentTime = get<0>(RisePileUp[i]);
			double currentEnergy = get<1>(RisePileUp[i]);
			int pileUpFlag = get<2>(RisePileUp[i]);

			while (i < RisePileUp.size() - 1)
			{

				if ((get<0>(RisePileUp[i + 1]) - currentTime) < 10000)
				{

					currentTime = (currentTime + get<0>(RisePileUp[i + 1])) / 2;
					currentEnergy = 0;
					pileUpFlag += get<2>(StoredEvents[i + 1]);
				}
				else
				{

					break;
				}

				i++;
			}

			TailPileUp.push_back(make_tuple(currentTime, currentEnergy, pileUpFlag));
		}

		cout << endl;
		cout << "Final list " << endl;

		for (unsigned int i = 0; i < TailPileUp.size(); i++)
		{

			cout << get<0>(TailPileUp[i]) << " " << get<1>(TailPileUp[i]) << " " << get<2>(TailPileUp[i]) << " " << G4Threading::G4GetThreadId() << endl;
		}

		sleep(3);
		StoredEvents = TailPileUp;
		return;
	}
}

void AppsDigitizer::ClearStoredEvents()
{

	StoredEvents.clear();
}
