define(['jquery', 'app/ui/layout', 'app/ui/detail', 'bootstrap'], function($, Layout, Detail) {

    var layout, detail, $search,
        $toggleProfile, $profileDropDown,
        $toggleCluster, $clusterSelection, selectedCluster,
        $toggleTimeSpan, $timeSpanSelection;

    function toggleClass($element, cssClass) {
        if ($element.hasClass(cssClass)) {
            $element.removeClass(cssClass);
            return false;
        }
        $element.addClass(cssClass);
        return true;
    }

    $(function() {

        layout = new Layout({debug: false});
        layout.init();

        detail = new Detail();

        $search = $('#search');

        $toggleProfile = $('#toggle-profile');
        $profileDropDown = $('#profile-dropdown');

        $toggleCluster = $('#toggle-cluster');
        $clusterSelection = $('#cluster-selection');

        $toggleTimeSpan = $('#toggle-timespan');
        $timeSpanSelection = $('#timespan-selection');

        // Profile
        $toggleProfile.on('click', function(e) {
            if (toggleClass($(e.delegateTarget), 'active')) {
                $profileDropDown.show();
            } else {
                $profileDropDown.hide();
            }
        });
        $profileDropDown.on('click', function(e) {
            e.stopPropagation();
        });

        // Toggle pinned
        $('#toggle-pinned').on('click', function(e) {
            toggleClass($(e.delegateTarget), 'active');
        });

        // Toggle search bar
        $('#toggle-search, #close-search').on('click', function() {
            toggleClass($search, 'opened');
        });

        // Search bar active
        $('#search-input')
            .on('focus', function(e) {
                e.stopPropagation();
                $(e.delegateTarget).parents('form').addClass('active');
            })
            .on('blur', function(e) {
                e.stopPropagation();
                $(e.delegateTarget).parents('form').removeClass('active');
            });

        // Toggle selector
        $('.selector button').on('click', function(e) {
            //e.stopPropagation();
            var $selector = $(e.target).parents('.selector').eq(0);
            toggleClass($selector, 'opened');
        });

        // Cluster selection
        $('#cluster').on('click', '.choices .cluster-value', function(e) {
            //e.stopPropagation();
            var $target = $(e.target);
            if (!$target.hasClass('cluster-value')) {
                $target = $target.parents('.cluster-value').eq(0);
            }

            var selection = $target.data('cluster');
            $clusterSelection.removeClass('cluster-' + selectedCluster);
            selectedCluster = selection;
            if (selectedCluster == 'none') {
                $clusterSelection.hide();
                $toggleCluster.find('.cluster-icon').show();
            } else {
                $clusterSelection.css({display: 'inline-block'}).addClass('cluster-' + selectedCluster);
                $toggleCluster.find('.cluster-icon').hide();
            }

            $toggleCluster.trigger('click');
        });

        // Timespan selection
        $('#timespan').on('click', '.choices a', function(e) {
            //e.stopPropagation();
            var value = (function(selection) {
                switch (selection) {
                    case 0 : return 'L';
                    case 1 : return 'W';
                    case 2 : return '1m';
                    case 3 : return '2m';
                }
            })($(e.target).closest('a').data('timespan'));

            $timeSpanSelection.html(value);
            $toggleTimeSpan.trigger('click');
        });

        // Filters
        $('.filter-toggle').on('click', function(e) {
            var $group = $(e.delegateTarget).parents('.filter-group').eq(0);
            if (toggleClass($group, 'active')) {
                $group.find('.filter-filters').collapse('show');
            } else {
                $group.find('.filter-filters').collapse('hide');
            }
        });
        $('.filter-group ul a').on('click', function(e) {
            toggleClass($(e.delegateTarget), 'active');
        });

        $('.filter-group.active .collapse').collapse('show');

        $('.filter-group .collapse')
            .on('hidden.bs.collapse', function (e) {
                e.stopPropagation();
                // TODO filtersScroll.redraw();
            })
            .on('shown.bs.collapse', function (e) {
                e.stopPropagation();
                // TODO filtersScroll.redraw();
            });

        // Thumbs
        $('.document, #detail')
            .on('click', '.thumb-pin', function(e) {
                e.stopPropagation();
                toggleClass($(e.target).closest('.thumb-pin'), 'active');
            })
            .on('click', '.thumb-remove', function(e) {
                e.stopPropagation();
                // TODO
            });

        // Detail
        $('.document')
            .on('click', '.thumb .title a', function(e) {
                e.preventDefault();
                detail.load();
                return false;
            });
        $('#detail-close, #backdrop').on('click', function() {
            detail.close();
        });
        $('#detail-next, #detail-prev').on('click', function() {
            detail.load();
        });

        $(window)
            // Close opened elements on window click
            .on('click', function(e) {
                e.stopPropagation();

                // Profile dropdown
                if (0 == $(e.target).closest('#toggle-profile').length) {
                    $('#toggle-profile').removeClass('active');
                    $('#profile-dropdown').hide();
                }

                // Selectors
                if (0 == $(e.target).closest('.selector').length) {
                    $('.selector').removeClass('opened');
                }
            })

    });
});

