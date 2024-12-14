from particle import*

def init_sph(space):
    for _ in range(100):
            x = np.random.uniform(space.bounds.left, space.bounds.right)
            y = np.random.uniform(space.bounds.top, space.bounds.bottom)
            space.particles.append(Particle(x, y))

def main():
    pygame.init()
    screen = pygame.display.set_mode((800, 600))
    space = Fluid_Space((200,200),(200,200))
    clock = pygame.time.Clock()
    flame = Flame(300, 450)
    
    init_sph(space)
    start_time = pygame.time.get_ticks()

    running = True
    while running:
        screen.fill((0, 0, 0))
        dt = clock.tick(60) / 1000
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

        elapsed_time = (pygame.time.get_ticks() - start_time) / 1000
        space.update()
        pygame.draw.rect(screen,'white',space.bounds,2)
        space.draw_particles(screen)
        if flame!=None:
            flame.apply_heat_to_particles(space, dt)
            flame.draw_flame(screen,dt)
        if elapsed_time > 15 :
            flame = None
        pygame.display.flip()
        clock.tick(60)

    pygame.quit()

if __name__ == "__main__":
    main()
