# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

function assert_particles_in_mesh(particles :: ParticleData, mesh :: Mesh)
    for (part_x, cell_i) in zip(particles.pos, particles.cell)
        if !cell_contains(mesh.cells[cell_i], part_x)
            error("Particle at $part_x is not inside of cell $cell_i ($(mesh.cells[cell_i]))")
        end
    end
end


function assert_cell_part_count(particles :: ParticleData, cell_part_count :: Vector)
    actual_cell_part_count = zeros(length(cell_part_count))
    for cell in particles.cell
        actual_cell_part_count[cell] += 1
    end
    for i in 1:length(cell_part_count)
        if actual_cell_part_count[i] != cell_part_count[i]
            error("Particle count in cell $i is inconsistent n=$(cell_part_count[i]), expected n=$(actual_cell_part_count[i])")
        end
    end
end