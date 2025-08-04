// multiple rider case attempt
/*
Passenger List: elevator maintains list of passengers (pairs of person_id and destination floor)
Batch Processing: elevator picks up as many passengers as possible from the queue if they are on the same floor, respecting the maximum occupancy
Route Management: elevator manages passengers' routes, moving to each passenger's destination floor and logging the necessary information
Complete Logging: detailed logs for each step including entering and exiting the elevator are added to provide comprehensive traceability
also
Occupancy Declaration: occupancy variable initialized to 0 at the start of the elevator function
Occupancy Management: increment occupancy each time a passenger enters the elevator and decrement each time a passenger exits.
*/

#ifndef HW6_ELEVATOR_HPP
#define HW6_ELEVATOR_HPP

#include <iostream>
#include <thread>
#include <mutex>
#include <queue>
#include <chrono>
#include <random>
#include <atomic>
#include <vector>
#include <condition_variable>
#include <tuple>
#include <unordered_map>
#include <algorithm>

using namespace std;

// Constants
const int NUM_FLOORS = 50;           // Total floors in the building
const int NUM_ELEVATORS = 6;         // Number of elevators
const int MAX_OCCUPANCY = 5;         // Maximum passengers per elevator
const int MAX_WAIT_TIME = 5000;      // Max wait time for person generation (ms)

// Global synchronization primitives
mutex cout_mtx;                      // Mutex for thread-safe cout
mutex queue_mtx;                     // Mutex for global queue (unused in current logic)
condition_variable cv;               // Condition variable for elevator wake-up

// Global data structures
queue<tuple<int, int, int>> global_queue;        // Global request queue (retained but unused)
vector<queue<tuple<int, int, int>>> floor_queues(NUM_FLOORS);  // Per-floor request queues
mutex floor_mtx[NUM_FLOORS];                     // Mutexes for each floor's queue
vector<int> elevator_positions(NUM_ELEVATORS, 0); // Current floor of each elevator
atomic<int> num_people_serviced(0);              // Total completed passenger requests
vector<int> global_passengers_serviced(NUM_ELEVATORS, 0); // Passengers served per elevator
int npeople;                                     // Total number of people to service
atomic<bool> all_done(false);

// User-added globals
struct Request {
    int person_id;
    int from;
    int to;
    bool assigned;
    Request(int pid, int f, int t) : person_id(pid), from(f), to(t), assigned(false) {}
};

struct Passenger {
    int person_id;
    int destination;
};

struct ElevatorState {
    int current_floor = 0;
    int direction = 0;
    vector<Passenger> passengers;
    int target_floor = -1;
    bool active = false;
};

static vector<Request> request_queue;
static mutex request_mutex;
static condition_variable request_cv;
static ElevatorState elevator_states[NUM_ELEVATORS];

void person(int id) {
    default_random_engine gen(random_device{}());
    uniform_int_distribution<int> dist(0, NUM_FLOORS - 1);
    int from = dist(gen);
    int to = dist(gen);
    while (to == from) to = dist(gen);

    {
        lock_guard<mutex> lock(request_mutex);
        request_queue.emplace_back(id, from, to);
    }

    {
        lock_guard<mutex> lock(cout_mtx);
        cout << "Person " << id << " wants to go from floor " << from << " to floor " << to << endl << flush;
    }
    request_cv.notify_all();
}

void elevator(int id) {
    ElevatorState& state = elevator_states[id];
    int& serviced = global_passengers_serviced[id];

    while (true) {
        if (all_done) {
            lock_guard<mutex> lock(cout_mtx);
            cout << "Elevator " << id << " has finished servicing all people." << endl;
            cout << "Elevator " << id << " serviced "
                 << global_passengers_serviced[id] << " passengers." << endl;
            return;
        }

        unique_lock<mutex> lock(request_mutex);
        request_cv.wait(lock, [] {
            return !request_queue.empty() || all_done;
        });

        if (all_done) {
            lock_guard<mutex> lock(cout_mtx);
            cout << "Elevator " << id << " has finished servicing all people." << endl;
            cout << "Elevator " << id << " serviced "
                 << global_passengers_serviced[id] << " passengers." << endl;
            return;
        }

        if (!state.active || state.passengers.empty()) {
            int min_dist = NUM_FLOORS + 1;
            int chosen_idx = -1;
            for (int i = 0; i < request_queue.size(); ++i) {
                Request& r = request_queue[i];
                if (!r.assigned) {
                    int dist = abs(r.from - state.current_floor);
                    if (dist < min_dist) {
                        min_dist = dist;
                        chosen_idx = i;
                    }
                }
            }

            if (chosen_idx != -1) {
                Request& r = request_queue[chosen_idx];
                r.assigned = true;
                state.active = true;
                state.target_floor = r.from;
                state.direction = (r.from > state.current_floor) ? 1 : ((r.from < state.current_floor) ? -1 : 0);
            } else {
                continue;
            }
        }
        lock.unlock();

        while (state.active) {
            this_thread::sleep_for(chrono::milliseconds(50));
            int next = state.current_floor + state.direction;
            if (next < 0 || next >= NUM_FLOORS) {
                state.direction = -state.direction;
                continue;
            }

            {
                lock_guard<mutex> lock(cout_mtx);
                cout << "Elevator " << id << " moving from floor " << state.current_floor
                     << " to floor " << next << endl << flush;
            }
            state.current_floor = next;

            // Drop off passengers
            auto drop_it = remove_if(state.passengers.begin(), state.passengers.end(), [&](Passenger& p) {
                if (p.destination == state.current_floor) {
                    lock_guard<mutex> lock(cout_mtx);
                    cout << "Person " << p.person_id << " arrived at floor "
                         << state.current_floor << endl << flush;
                    ++serviced;
                    ++num_people_serviced;
                    return true;
                }
                return false;
            });
            state.passengers.erase(drop_it, state.passengers.end());

            // Pick up passengers
            {
                lock_guard<mutex> lock2(request_mutex);
                for (auto it = request_queue.begin(); it != request_queue.end();) {
                    Request& r = *it;
                    bool can_pickup = false;

                    if (r.assigned && r.from == state.current_floor) {
                        if (state.passengers.empty()) {
                            can_pickup = true;
                        } else if (state.direction == 0) {
                            can_pickup = true;
                        } else if (state.direction == 1 && r.to > r.from) {
                            can_pickup = true;
                        } else if (state.direction == -1 && r.to < r.from) {
                            can_pickup = true;
                        }
                    }

                    if (can_pickup && (int)state.passengers.size() < MAX_OCCUPANCY) {
                        state.passengers.push_back({r.person_id, r.to});
                        {
                            lock_guard<mutex> lock(cout_mtx);
                            cout << "Person " << r.person_id << " entered elevator " << id << endl << flush;
                        }
                        if (state.direction == 0) {
                            state.direction = (r.to > state.current_floor) ? 1 : -1;
                        }
                        if (state.direction == 1 && r.to > state.target_floor) state.target_floor = r.to;
                        if (state.direction == -1 && r.to < state.target_floor) state.target_floor = r.to;
                        it = request_queue.erase(it);
                        continue;
                    }

                    ++it;
                }
            }

            // Update target
            if (!state.passengers.empty()) {
                if (state.direction == 1) {
                    int max_dest = state.target_floor;
                    for (auto& p : state.passengers)
                        if (p.destination > max_dest) max_dest = p.destination;
                    state.target_floor = max_dest;
                } else if (state.direction == -1) {
                    int min_dest = state.target_floor;
                    for (auto& p : state.passengers)
                        if (p.destination < min_dest) min_dest = p.destination;
                    state.target_floor = min_dest;
                }
            }

            if (state.passengers.empty() && state.current_floor == state.target_floor) {
                state.active = false;
                state.direction = 0;

                if (num_people_serviced >= npeople) {
                    all_done = true;
                    request_cv.notify_all(); // wake up other elevators
                }

                break;
            }
        }
    }
}

#endif // ELEVATOR_HPP