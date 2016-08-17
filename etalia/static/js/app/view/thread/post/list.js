define([
    'app',
    'text!app/templates/thread/post/list.hbs',
    'app/view/thread/post/thumb',
    'app/view/thread/post/form'
], function (App, template) {

    App.View.Thread.PostList = App.Backbone.View.extend({

        tagName: 'div',
        className: 'thread-posts',

        template: App.Handlebars.compile(template),

        thread: null,
        form: null,

        initialize: function (options) {
            if (!options.thread) {
                throw 'Expected instance of App.Model.Thread';
            }

            this.thread = options.thread;

            this.listenTo(this.model, "sync", this.render);
        },

        renderForm: function () {
            if (this.form) {
                this.form.remove();
            }

            this.form = App.View.Thread.PostForm.create({}, {
                $target: this.$('[data-post-add-form]')
            });

            this.listenToOnce(this.form, 'validation_success', this.submitForm);
        },

        submitForm: function () {
            var that = this;
            this.form.model.save({
                user: App.getCurrentUser(),
                thread: this.thread
            }, {
                success: function () {
                    that.renderForm();
                },
                error: function() {
                    that.form.model.destroy();
                }
            });
        },

        render: function () {
            App.log('PostListView::render');

            var is_member = this.thread.isMember(App.getCurrentUser());

            // TODO memberOfThread
            this.$el.html(this.template({
                is_member: is_member
            }));

            var that = this,
                $list = this.$('.thread-posts-list');
            this.model.each(function (post) {
                that.pushSubView(
                    App.View.Thread.PostThumb.create({
                        model: post
                    }, {
                        $target: $list,
                        append: true
                    })
                );
            });

            this.pushSubView(
                App.View.User.Thumb.create({
                    model: App.getCurrentUser()
                }, {
                    $target: this.$('.user-placeholder')
                })
            );

            if (is_member) {
                // TODO handle subView
                this.renderForm();
            }

            return this;
        }
    });

    App.View.Thread.PostList.create = function (options, createOptions) {
        options = options || {};
        if (!options.thread) {
            throw 'options.thread is expected to be instance of App.Model.Thread';
        }
        options.model = options.thread.get('posts');

        var list = new App.View.Thread.PostList(options);

        if (createOptions) {
            App.View.create(list, createOptions);
        }

        return list;
    };

    return App.View.Thread.PostList;
});
