// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

/******** Functions to find iProc starts ********/

// --------------------------------------------------------------------------
// Compute iProc on a plane divided first into left and right, then into up and down
// i.e. 2 --- 3
//      |     |
//      0 --- 1
// See quadtree.cpp:133-134 in function Quadtree::new_node for details
// --------------------------------------------------------------------------

int pos_iProc_surface(precision_t oLR,
                      precision_t oDU,
                      precision_t dLR,
                      precision_t dDU,
                      precision_t LR,
                      precision_t DU,
                      int numProcs) {
    // The iProc of the lower left corner
    int retProc = 0;

    // Repeat until there is only one processor in the described area
    while (numProcs != 1) {
        // Divide into 4 parts
        dLR /= 2;
        dDU /= 2;

        if (LR >= oLR + dLR) {
            // The point is on the right, and the iProc of the lower left
            // corner where the point is located is [1/4] or [3/4]
            // * remainProc, depending on DU
            retProc += numProcs / 4;
            // Update the origin
            oLR += dLR;
        }

        if (DU >= oDU + dDU) {
            // The point is above
            retProc += numProcs / 2;
            // Update the origin
            oDU += dDU;
        }

        numProcs /= 4;
    }

    return retProc;
}

// --------------------------------------------------------------------------
// Compute iProc in sphere for the given point
// --------------------------------------------------------------------------

int pos_iProc_sphere(precision_t lon_in, precision_t lat_in) {
    return pos_iProc_surface(Sphere::ORIGINS(0) * cPI,
                             Sphere::ORIGINS(1) * cPI,
                             Sphere::RIGHTS(0) * cPI,
                             Sphere::UPS(1) * cPI,
                             lon_in,
                             lat_in,
                             nGrids);
}

// --------------------------------------------------------------------------
// Return true if c is in [a, a + b) or (a + b, a]
// --------------------------------------------------------------------------

bool inRange(precision_t a, precision_t b, precision_t c) {
    if (b > 0) {
        return c >= a && c < a + b;
    } else {
        return c <= a && c > a + b;
    }
}

// --------------------------------------------------------------------------
// Compute iProc in cubesphere for the given point
// --------------------------------------------------------------------------

int pos_iProc_cube(precision_t lon_in, precision_t lat_in) {
    // Transfer polar coordinate to cartesian coordinate
    arma_vec point = sphere_to_cube(lon_in, lat_in);

    // Store the index of surface, LR and DU
    // INITIAL VALUES ARE HARD-CODED BASED ON SURFACE 5(6) TO
    // SOLVE BOUNDARY PROBLEMS
    int iFace = 0, iLR = 1, iDU = 0;
    bool on_surface = false;
    // Find the surface where the point is located
    for (;iFace < 6 && !on_surface; ++iFace) {
        // Assume that the point is on current surface
        on_surface = true;
        int sLR = -1, sDU = -1;

        // Check three dimensions
        for (int j = 0; j < 3 && on_surface; ++j) {
            if (CubeSphere::RIGHTS(iFace, j) != 0) {
                sLR = j;
                on_surface = inRange(CubeSphere::ORIGINS(iFace, j),
                                     CubeSphere::RIGHTS(iFace, j),
                                     point(j));
            } else if (CubeSphere::UPS(iFace, j) != 0) {
                sDU = j;
                on_surface = inRange(CubeSphere::ORIGINS(iFace, j),
                                     CubeSphere::UPS(iFace, j),
                                     point(j));
            } else {
                on_surface = (CubeSphere::ORIGINS(iFace, j) == point(j));
            }
        }

        // Assign iLR and iDU if the point is indeed on current surface
        if (on_surface) {
            iLR = sLR;
            iDU = sDU;
        }
    }
    // Minus 1 due to for loop mechanism
    --iFace;

    // Change the sign if delta is negative
    precision_t coefLR, coefDU;
    coefLR = CubeSphere::RIGHTS(iFace, iLR) > 0 ? 1.0 : -1.0;
    coefDU = CubeSphere::UPS(iFace, iDU) > 0 ? 1.0 : -1.0;
    // Compute iProc in the surface
    int ip = pos_iProc_surface(CubeSphere::ORIGINS(iFace, iLR) * coefLR,
                               CubeSphere::ORIGINS(iFace, iDU) * coefDU,
                               CubeSphere::RIGHTS(iFace, iLR) * coefLR,
                               CubeSphere::UPS(iFace, iDU) * coefDU,
                               point(iLR) * coefLR,
                               point(iDU) * coefDU,
                               nGrids / 6);
    // Consider surface number
    return ip + iFace * nGrids / 6;
}

/******** Functions to find iproc ends ********/


/******** Functions related to MPI starts ********/

// --------------------------------------------------------------------------
// MPI message structure
// --------------------------------------------------------------------------

template<typename T>
struct mpi_msg_t {
    // The message to send
    std::vector<T> buf;
    // The processor to send
    int rank;
    // Default constructor
    mpi_msg_t()
        : buf(), rank() {}
    // Move constructor
    mpi_msg_t(std::vector<T> &&_buf, const int _rank)
        : buf(std::move(_buf)), rank(_rank) {}
};

// --------------------------------------------------------------------------
// Debug
// --------------------------------------------------------------------------

template<typename T>
std::ostream& operator<<(std::ostream &out,
                         const std::vector<mpi_msg_t<T>> &msg) {
    for (auto& it : msg) {
        out << "\tTo " << it.rank << ":\n";
        for (auto& it2 : it.buf) {
            out << "\t\t" << it2 << '\n';
        }
        out << '\n';
    }
    out << '\n';
    return out;
}

/**
 * \brief Send the messages to their destination, receive
 *        messages from other processors
 *        The main problem solved is each processor only
 *        hold the messages to send but doesn't know how
 *        many messages will come from which processor
 * \param msg_in The messages to send
 * \param msg_out The messages received
 * \param comm MPI communicator
 * \return True if the function succeeds
 */
template<typename T>
bool CUSTOM_MPI_SENDRECV(const std::vector<mpi_msg_t<T>> &msg_in,
                         std::vector<mpi_msg_t<T>> &msg_out,
                         const MPI_Comm &comm) {
    // Promise previous send and receive have all ended
    MPI_Barrier(comm);

    // Clear the output variables
    msg_out.clear();

    // Get the rank the size of communicator
    // Processor with rank 0 checks whether all ends
    const int num = msg_in.size();
    int comm_rank, comm_size;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);

    // Define constants
    const int TAG_RESERVED = 3;
    const char ALL_FINISH_TAG = 0;
    const char PROC_FINISH_TAG = 1;
    const char ACK_RECV_TAG = 2;

    // Only used to check correctness
    std::vector<MPI_Request> request;

    // Send
    for (int i = 0; i < num; ++i) {
        request.emplace_back();
        MPI_Isend(msg_in[i].buf.data(),
                  msg_in[i].buf.size() * sizeof(T),
                  MPI_BYTE,
                  msg_in[i].rank,
                  i + TAG_RESERVED,
                  comm,
                  &request.back());
    }

    // Receive
    // Record how many messages are sent successfully and
    // how many processors finish sending messages
    // finished_proc is only used by iProc 0
    int success_sent = 0, finished_proc = 0;
    bool done = false;
    // If nothing needs to send, then the processor already finish sending
    if (!num) {
        if (!comm_rank) {
            // Increase finished_proc if iProc = 0
            ++finished_proc;
        } else {
            // Send signal to iProc 0 otherwise
            request.emplace_back();
            MPI_Isend(NULL, 0, MPI_BYTE, 0, PROC_FINISH_TAG,
                      comm, &request.back());
        }
    }
    while (!done) {
        MPI_Status status;
        int count;
        // Probe the size of next message
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
        MPI_Get_count(&status, MPI_BYTE, &count);
        if (count) {
            // Receiving actual messages
            mpi_msg_t<T> msg;
            // Set rank and count, create space to receive message
            msg.buf.resize(count / sizeof(T));
            msg.rank = status.MPI_SOURCE;
            MPI_Recv(msg.buf.data(),
                     count,
                     MPI_BYTE,
                     status.MPI_SOURCE,
                     status.MPI_TAG,
                     comm,
                     MPI_STATUS_IGNORE);
            msg_out.push_back(msg);
            // Send ack message to the origin of this message
            request.emplace_back();
            MPI_Isend(NULL, 0, MPI_BYTE, status.MPI_SOURCE,
                      ACK_RECV_TAG, comm, &request.back());
        } else {
            // The message is a control signal with size = 0
            MPI_Recv(NULL, 0, MPI_BYTE, status.MPI_SOURCE,
                     status.MPI_TAG, comm, &status);
            switch (status.MPI_TAG) {
                case ACK_RECV_TAG: {
                    // Increase the success_sent
                    ++success_sent;
                    if (success_sent == num) {
                        // The processor finish sending
                        if (!comm_rank) {
                            // Increase finished_proc for rank 0
                            ++finished_proc;
                        } else {
                            // Send proc finish to rank 0 for other processors
                            request.emplace_back();
                            MPI_Isend(NULL, 0, MPI_BYTE, 0, PROC_FINISH_TAG,
                                      comm, &request.back());
                        }
                    }
                    break;
                }
                case PROC_FINISH_TAG: {
                    // Only processor with rank 0 will get this tag
                    ++finished_proc;
                    break;
                }
                case ALL_FINISH_TAG: {
                    // Exit the while loop
                    done = true;
                    break;
                }
                default: {
                    // This should never be executed
                    return false;
                }
            }
            // Check whether everything is done by processor with rank 0
            if (finished_proc == comm_size) {
                // Send all finish to all other processors
                for (int i = 1; i < comm_size; ++i) {
                    request.emplace_back();
                    MPI_Isend(NULL, 0, MPI_BYTE, i, ALL_FINISH_TAG,
                              comm, &request.back());
                }
                done = true;
            }
        }
    }

    // Correctness check using request
    for (auto& it : request) {
        int flag;
        MPI_Test(&it, &flag, MPI_STATUS_IGNORE);
        if (!flag) {
            // This should never be executed
            return false;
        }
    }

    // Promise all send and receive end
    MPI_Barrier(comm);

    return true;
}

/**
 * \brief Send and receive messages. The origin and number
 *        of messages are determined
 *        BY DESIGN THE DESTINATION OF ANY TWO MESSAGES IN
 *        msg_in CANNOT BE EQUAL
 * \param msg_in The messages to send
 * \param msg_out The messages received
 * \param comm MPI communicator
 * \return True if the function succeeds
 */
template<typename T>
bool CUSTOM_MPI_SENDRECV2(const std::vector<mpi_msg_t<T>> &msg_in,
                          std::vector<mpi_msg_t<T>> &msg_out,
                          const MPI_Comm &comm) {
    // Promise previous send and receive have all ended
    MPI_Barrier(comm);

    // Only used to check correctness
    std::vector<MPI_Request> request;

    // Send
    for (auto& it : msg_in) {
        request.emplace_back();
        MPI_Isend(it.buf.data(),
                  it.buf.size() * sizeof(T),
                  MPI_BYTE,
                  it.rank,
                  0,
                  comm,
                  &request.back());
    }

    // Receive
    for (auto& it : msg_out) {
        MPI_Recv(it.buf.data(),
                 it.buf.size() * sizeof(T),
                 MPI_BYTE,
                 it.rank,
                 0,
                 comm,
                 MPI_STATUS_IGNORE);
    }

    // Promise all send and receive end
    MPI_Barrier(comm);

    // Correctness check using request
    for (auto& it : request) {
        int flag;
        MPI_Test(&it, &flag, MPI_STATUS_IGNORE);
        if (!flag) {
            // This should never be executed
            return false;
        }
    }

    return true;
}

/******** Functions related to MPI ends ********/

// --------------------------------------------------------------------------
// Structure of longitude and latitude
// --------------------------------------------------------------------------

struct lon_lat_t {
    precision_t lon;
    precision_t lat;
};

// --------------------------------------------------------------------------
// Debug
// --------------------------------------------------------------------------

std::ostream& operator<<(std::ostream &out, const lon_lat_t &lonlat) {
    out << "lon = " << lonlat.lon << "\tlat = " << lonlat.lat << '\n';
    return out;
}

// --------------------------------------------------------------------------
// Wrap the point only for spherical grid
// Return -1 if the wrap crosses the north or south pole, 1 otherwise
// --------------------------------------------------------------------------

int wrap_point_sphere(precision_t &lon, precision_t &lat) {
    int inverse = 1;
    if (lat < -0.5 * cPI) {
        lat = -cPI - lat;
        inverse = -1;
        lon += cPI;
    } else if (lat > 0.5 * cPI) {
        lat = cPI - lat;
        inverse = -1;
        lon += cPI;
    }
    while (lon < 0) {
        lon += cTWOPI;
    }
    while (lon > cTWOPI) {
        lon -= cTWOPI;
    }
    return inverse;
}

// --------------------------------------------------------------------------
// Initialize the connection between processors for the geo grid
// --------------------------------------------------------------------------

void Grid::init_connection() {
    // Initialize communicator. Message exchange is limited in each
    // ensemble and rank 0 in each ensemble has special purpose
    MPI_Comm_split(aether_comm, iMember, iGrid, &grid_comm);

    // Classify the {lon, lat} for all ghost cells based on destination
    std::unordered_map<int, std::vector<lon_lat_t>> lonlatmap;
    for (int64_t iLon = 0; iLon < nLons; ++iLon) {
        for (int64_t iLat = 0; iLat < nLats; ++iLat) {
            // Read the lon, lat of ghost cell
            precision_t lon = geoLon_scgc(iLon, iLat, 0);
            precision_t lat = geoLat_scgc(iLon, iLat, 0);
            // Find which processor can handle it and whether
            // the message crosses the north/south pole
            int dest, inverse = 1;
            if (IsCubeSphereGrid) {
                dest = pos_iProc_cube(lon, lat);
            } else {
                // Spherical grid needs wrap
                inverse = wrap_point_sphere(lon, lat);
                dest = pos_iProc_sphere(lon, lat);
            }

            // Prepare to send the lon, lat to dest processor
            lonlatmap[dest].push_back({lon, lat});
            // Record where to place data when receiving from that processor
            exch_recv[dest].push_back({iLon, iLat, inverse});
            // Skip non ghost cells
            if (iLat == nGCs - 1 && iLon >= nGCs && iLon < nLons - nGCs) {
                iLat = nLats - nGCs - 1;
            }
        }
    }

    // Construct send messages
    std::vector<mpi_msg_t<lon_lat_t>> msg_in;
    for (auto& it : lonlatmap) {
        msg_in.emplace_back(std::move(it.second), it.first);
    }
    // Create space to receive messages
    std::vector<mpi_msg_t<lon_lat_t>> msg_out;

    // Send all longitude and latitude to where it can be handled
    // Receive all longitude and latitude this processor needs to handle
    CUSTOM_MPI_SENDRECV(msg_in, msg_out, grid_comm);
    
    // Build interpolation coefficients based on lon and lat
    for (auto& it : msg_out) {
        // Initialize the inputs to set_interp_coefs
        std::vector<precision_t> Lons;
        std::vector<precision_t> Lats;
        std::vector<precision_t> Alts;

        // Combine each altitude with lat_lon into points
        for (int64_t iAlt = nGCs; iAlt < nAlts - nGCs; ++iAlt) {
            for (auto& it2 : it.buf) {
                Lons.push_back(it2.lon);
                Lats.push_back(it2.lat);
                Alts.push_back(geoAlt_scgc(0, 0, iAlt));
            }
        }
        // Compute interpolation coefficients
        set_interpolation_coefs(Lons, Lats, Alts);

        // Fix precision issues
        for (auto& it : interp_coefs) {
            if (!it.in_grid) {
                exit(1);
            }
            if (it.rRow < 0.0001) {
                it.rRow = 0;
            } else if (it.rRow > 0.9999) {
                it.rRow = 1;
            }
            if (it.rCol < 0.0001) {
                it.rCol = 0;
            } else if (it.rCol > 0.9999) {
                it.rCol = 1;
            }
        }
        // Store the coefficients
        exch_send[it.rank] = std::move(interp_coefs);
    }
}

// --------------------------------------------------------------------------
// Update ghost cells with values from other processors
// --------------------------------------------------------------------------

void Grid::exchange(arma_cube &data, const bool pole_inverse) {
    // Construct send and receive buffer
    std::vector<mpi_msg_t<precision_t>> msg_send;
    std::vector<mpi_msg_t<precision_t>> msg_recv;

    // Do the interpolation for other processors
    for (auto& it : exch_send) {
        std::swap(interp_coefs, it.second);
        msg_send.emplace_back(get_interpolation_values(data),
                              it.first);
        std::swap(interp_coefs, it.second);
    }

    // Initialize the size of receive buffer
    const int64_t numalt = nAlts - 2 * nGCs;
    for (auto& it : exch_recv) {
        msg_recv.emplace_back(std::vector<precision_t>(it.second.size()*numalt),
                              it.first);
    }

    // Exchange message
    CUSTOM_MPI_SENDRECV2(msg_send, msg_recv, grid_comm);

    // Fill the ghost cell with received message
    for (auto& it : msg_recv) {
        const std::vector<idx2d_t> &idx = exch_recv[it.rank];
        for (int64_t iAlt = nGCs; iAlt < nAlts - nGCs; ++iAlt) {
            for (size_t i = 0; i < idx.size(); ++i) {
                if (pole_inverse) {
                    it.buf[(iAlt - nGCs) * idx.size() + i] *= idx[i].inverse;
                }
                data(idx[i].ilon, idx[i].ilat, iAlt) =
                    it.buf[(iAlt - nGCs) * idx.size() + i];
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Exchange messages between processors
// -----------------------------------------------------------------------------

bool Neutrals::exchange(Grid &grid) {
    // For each species, exchange if its DoAdvect is true
    for (int i = 0; i < nSpecies; ++i) {
        if (species[i].DoAdvect) {
            grid.exchange(species[i].density_scgc, false);
        }
    }
    // Exchange temperature
    grid.exchange(temperature_scgc, false);
    // Exchange velocity
    grid.exchange(velocity_vcgc[0], true);
    grid.exchange(velocity_vcgc[1], true);
    grid.exchange(velocity_vcgc[2], false);
    return true;
}
