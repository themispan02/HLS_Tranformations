/*
 * Gravitational N-body simulation (float): plain C reference without HLS
 * directives. Exploration, not used in the thesis results; the thesis uses the
 * MachSuite molecular dynamics kernel instead.
 *
 * Each time step updates the positions, computes all-pairs gravity (O(N^2))
 * and updates the velocities (symplectic Euler). Parameters follow the AMD AI
 * Engine N-body tutorial.
 */

#include <math.h>

#define N_BODIES       160
#define N_DIMS         3
#define TIMESTEPS      4

#define G_CONST        1.0f
#define SOFTENING_SQ   1000.0f
#define DT             1.0f


void nbody_naive(
    float pos_x[N_BODIES],
    float pos_y[N_BODIES],
    float pos_z[N_BODIES],
    float vel_x[N_BODIES],
    float vel_y[N_BODIES],
    float vel_z[N_BODIES],
    float mass[N_BODIES])
{
    int i, j;

    /* Ping-pong buffers: even steps read pos/vel and write buf_*, odd steps the reverse */
    float buf_pos_x[N_BODIES];
    float buf_pos_y[N_BODIES];
    float buf_pos_z[N_BODIES];
    float buf_vel_x[N_BODIES];
    float buf_vel_y[N_BODIES];
    float buf_vel_z[N_BODIES];

    float *src_px, *src_py, *src_pz;
    float *src_vx, *src_vy, *src_vz;
    float *dst_px, *dst_py, *dst_pz;
    float *dst_vx, *dst_vy, *dst_vz;

    Loop_timestep: for (int t = 0; t < TIMESTEPS; t++) {

        if (t % 2 == 0) {
            src_px = pos_x;   src_py = pos_y;   src_pz = pos_z;
            src_vx = vel_x;   src_vy = vel_y;   src_vz = vel_z;
            dst_px = buf_pos_x; dst_py = buf_pos_y; dst_pz = buf_pos_z;
            dst_vx = buf_vel_x; dst_vy = buf_vel_y; dst_vz = buf_vel_z;
        } else {
            src_px = buf_pos_x; src_py = buf_pos_y; src_pz = buf_pos_z;
            src_vx = buf_vel_x; src_vy = buf_vel_y; src_vz = buf_vel_z;
            dst_px = pos_x;   dst_py = pos_y;   dst_pz = pos_z;
            dst_vx = vel_x;   dst_vy = vel_y;   dst_vz = vel_z;
        }
            dst_px[i] = src_px[i] + src_vx[i] * DT;
            dst_py[i] = src_py[i] + src_vy[i] * DT;
            dst_pz[i] = src_pz[i] + src_vz[i] * DT;
        }


        Loop_i: for (i = 0; i < N_BODIES; i++) {

            float acc_x = 0.0f;
            float acc_y = 0.0f;
            float acc_z = 0.0f;

            float pi_x = dst_px[i];
            float pi_y = dst_py[i];
            float pi_z = dst_pz[i];

            Loop_j: for (j = 0; j < N_BODIES; j++) {

                float dx = dst_px[j] - pi_x;
                float dy = dst_py[j] - pi_y;
                float dz = dst_pz[j] - pi_z;

                /* Softening keeps i == j and close pairs finite */
                float dist_sq = dx * dx + dy * dy + dz * dz + SOFTENING_SQ;

                /* 1 / ||r||³  =  (1/√dist_sq)³ */
                float inv_dist  = 1.0f / sqrtf(dist_sq);
                float inv_dist3 = inv_dist * inv_dist * inv_dist;

                /* acc_i += G * m_j * r_ij / ||r||³ */
                float s = G_CONST * mass[j] * inv_dist3;

                acc_x += s * dx;
                acc_y += s * dy;
                acc_z += s * dz;
            }

            /* v_{t+1} = v_t + acc * dt */
            dst_vx[i] = src_vx[i] + acc_x * DT;
            dst_vy[i] = src_vy[i] + acc_y * DT;
            dst_vz[i] = src_vz[i] + acc_z * DT;
        }

    } /* end Loop_timestep */

    /* After an odd number of steps the result is in buf_*: copy it back */
    if (TIMESTEPS % 2 != 0) {
        for (i = 0; i < N_BODIES; i++) {
            pos_x[i] = buf_pos_x[i];
            pos_y[i] = buf_pos_y[i];
            pos_z[i] = buf_pos_z[i];
            vel_x[i] = buf_vel_x[i];
            vel_y[i] = buf_vel_y[i];
            vel_z[i] = buf_vel_z[i];
        }
    }
}