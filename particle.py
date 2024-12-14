import pygame
import numpy as np

G = np.array([0.0, 10.0])  # Gravity
REST_DENS = 100.0           # Rest density
GAS_CONST = 100.0          # Gas constant
H = 20.0                    # Kernel radius
HSQ = H ** 2                # Kernel radius squared
MASS = 1                  # Particle mass
VISC = 100.0                # Viscosity constant
DT = 0.001                 # Time step
POLY6 = 4.0 / (np.pi * H ** 8)
SPIKY_GRAD = -45.0 / (np.pi * H ** 6)
VISC_LAP = 45.0 / (np.pi * H ** 5)


class Particle:
    def __init__(self, x, y):
        jitter = np.random.random()
        self.x = np.array([x+ jitter, y])  # Position
        self.v = np.zeros(2)               # Velocity
        self.f = np.zeros(2)               # Force
        self.rho = 0.0                     # Density
        self.p = 0.0                       # Pressure

        self.temperature = 20.0

        self.n_p = []
    def thermal_equilibrium(self, other):
        temp= self.temperature + other.temperature
        self.temperature = temp/2
        other.temperature = temp/2

class Fluid_Space:
    def __init__(self,  bounds_size, bound_pos):
        self.bounds_width, self.bounds_height = bounds_size

        self.bounds_x,self.bounds_y  = bound_pos
        self.bounds = pygame.Rect(self.bounds_x, self.bounds_y, self.bounds_width, self.bounds_height)
        self.particles = []
        self.space_temp = 20.0

    def create_grid(self, h):
        cell_size = h
        grid = {}
        for p in self.particles:
            # 셀 좌표 계산
            cell_x = p.x[0] // cell_size
            cell_y = p.x[1] // cell_size
            cell_key = (cell_x, cell_y)

            # 셀에 입자 추가
            if cell_key not in grid:
                grid[cell_key] = []
            grid[cell_key].append(p)

        return grid

    def find_near_particles(self,grid, h):
        cell_size = h
        for p in self.particles:
            cell_x = p.x[0] // cell_size
            cell_y = p.x[1] // cell_size

            p.n_p = []  
            for dx in range(-1, 2):  
                for dy in range(-1, 2):  
                    neighbor_cell = (cell_x + dx, cell_y + dy)
                    if neighbor_cell in grid:
                        for neighbor in grid[neighbor_cell]:
                                p.n_p.append(neighbor)

    def compute_density_pressure(self):
        for pi in self.particles:
            pi.rho = 0.0
            for pj in self.particles:
                rij = pj.x - pi.x
                r2 = np.dot(rij, rij)
                if r2 < HSQ:
                    pi.rho += MASS * POLY6 * (HSQ - r2) ** 3
            pi.p = GAS_CONST * (pi.rho - REST_DENS)
            
    def compute_forces(self):
        for pi in self.particles:
            fpress = np.zeros(2)
            fvisc = np.zeros(2)
            for pj in pi.n_p:
                rij = pj.x - pi.x
                r = np.linalg.norm(rij)
                if r < H and r>0:
                    # Pressure force
                    fpress += (-rij / r) * MASS * (pi.p + pj.p) / (2 * pj.rho) * SPIKY_GRAD * (H - r) ** 3
                    # Viscosity force
                    fvisc += VISC * MASS * (pj.v - pi.v) / pj.rho * VISC_LAP * (H - r)
            fgrav = G * MASS / pi.rho
            pi.f = fpress + fvisc + fgrav

    def integrate(self):
        for p in self.particles:
            p.v += DT * p.f / p.rho
            p.x += DT * p.v

            # Boundary conditions
            if p.x[0] - 4 < self.bounds.left:
                p.v[0] *= -0.5
                p.x[0] = self.bounds.left  + 4
            if p.x[0] + 4 > self.bounds.right:
                p.v[0] *= -0.5
                p.x[0] = self.bounds.right -4
            if p.x[1] -4 < self.bounds.top:
                p.v[1] *= -0.5
                p.x[1] = H + self. bounds.top +4
            if p.x[1]+4 > self.bounds.bottom:
                p.v[1] *= -0.5
                p.x[1] = self.bounds.bottom -4
    
    def compute_heat_transfer(self):
        for pi in self.particles:
            for pj in pi.n_p:
                if pi != pj: 
                    pi.thermal_equilibrium(pj)
            pi.temperature -= 0.5*(pi.temperature-self.space_temp)*DT


    def draw_particles(self, screen):
        for particle in self.particles:
            normalized_temp = (particle.temperature - 20) / (100 - 20)
            normalized_temp = max(0, min(normalized_temp, 1))
            color = (int(255 * normalized_temp),int(255 * (1 - normalized_temp)),50)
            pygame.draw.circle(screen, color, (int(particle.x[0]), int(particle.x[1])), 4)

    def update(self):
        grid = self.create_grid(6)
        self.find_near_particles(grid, 6)
        self.compute_density_pressure()
        self.compute_forces()
        self.compute_heat_transfer()
        self.integrate()

class FlameParticle:
    alpha_layer_qty = 2
    alpha_glow_const = 2

    def __init__(self, x, y, r=5, lifetime=1):
        self.x = x
        self.y = y
        self.r = r
        self.original_r = r
        self.lifetime = lifetime  
        self.alpha_layers = FlameParticle.alpha_layer_qty
        self.alpha_glow = FlameParticle.alpha_glow_const
        max_surf_size = 2 * self.r * self.alpha_layers * self.alpha_layers * self.alpha_glow
        self.surf = pygame.Surface((max_surf_size, max_surf_size), pygame.SRCALPHA)

    def update(self, dt):
        self.lifetime -= dt  
        self.y -= (7 - self.r) * dt * 10 
        self.x += np.random.uniform(-self.r, self.r) * dt  
        self.r = max(1, int(self.original_r * (self.lifetime / 2)))  

    def draw(self, screen):
        max_surf_size = 2 * self.r * self.alpha_layers * self.alpha_layers * self.alpha_glow
        self.surf = pygame.Surface((max_surf_size, max_surf_size), pygame.SRCALPHA)
        for i in range(self.alpha_layers, -1, -1):
            alpha = 255 - i * (255 // self.alpha_layers - 5)
            if alpha <= 0:
                alpha = 0
            radius = self.r * i * i * self.alpha_glow
            lifetime_ratio = max(0, self.lifetime) / 2.0 
            if lifetime_ratio > 0.7:  
                ratio = (lifetime_ratio - 0.7) / 0.7
                r = 255
                g = int(150 * ratio) 
                b = 0
            else:  
                ratio = lifetime_ratio / 0.7
                r = int(255 * ratio + (1 - ratio) * 50)
                g = int(150 * ratio + (1 - ratio) * 50)
                b = int((1 - ratio) * 50)

            color = (r, g, b, alpha)
            pygame.draw.circle(self.surf, color, (self.surf.get_width() // 2, self.surf.get_height() // 2), radius)
        screen.blit(self.surf, self.surf.get_rect(center=(self.x, self.y)))

class Flame:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.flame_intensity = 2
        self.flame_particles = []

    def apply_heat_to_particles(self, fluid_space, dt):
        for flame_particle in self.flame_particles:
            for particle in fluid_space.particles:
                distance = np.linalg.norm(np.array([flame_particle.x, flame_particle.y]) - np.array([particle.x[0], particle.x[1]]))
                if distance < 10: 
                    particle.temperature += 50 * dt 
                    particle.temperature = min(particle.temperature, 100.0)  

    def draw_flame(self,screen, dt):
        for i in range(1):
            self.flame_particles.append(
                FlameParticle(self.x + np.random.randint(-5, 5), self.y, np.random.randint(1, 5), lifetime=2)
            )
        for i in self.flame_particles:
            if i.lifetime <= 0:  
                self.flame_particles.remove(i)
                del i
                continue
            i.update(dt) 
            i.draw(screen)