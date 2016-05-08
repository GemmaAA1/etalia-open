define([
    'app',
    'text!app/templates/user/list.hbs',
    'app/view/user/thumb'
], function (App, template) {

    App.View.User = App.View.User || {};

    App.View.User.List = App.Backbone.View.extend({

        tagName: 'div',
        className: 'thread-members thread-box',

        template: App.Handlebars.compile(template),

        events: {
            //"click #thread-next-page": "onNextPageClick"
        },

        initialize: function () {
            //options = App.defaults(defaults, options);

            this.listenTo(this.model, "reset update", this.render);
        },

        render: function () {
            App.log('UserListView::render');

            this.$el.html(this.template({}));

            var that = this,
                $list = this.$('.thread-members-list');

            this.model.each(function(user) {
                that.pushSubView(
                    App.View.User.Thumb.create({
                        model: user
                    }, {
                        $target: $list,
                        append: true
                    })
                );
            });

            return this;
        }
    });

    App.View.User.List.create = function(options, createOptions) {
        options = options || {};
        if (!options.model) {
            options.model = new App.Collection.Users();
        }

        var list = new App.View.User.List(options);
        if (createOptions) {
            App.View.create(list, createOptions);
        }

        return list;
    };

    return App.View.User.List;
});
